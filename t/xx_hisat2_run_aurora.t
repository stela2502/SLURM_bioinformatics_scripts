#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 9;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my (
	$value,        @values,  $exp,      $max_jobs, @files,
	$bigwigTracks, @options, $coverage, $genome,   $outpath,
);

my $exec = $plugin_path . "/../bin/hisat2_run_aurora.pl";
ok( -f $exec, 'the script has been found' );
$outpath = "$plugin_path/data/output/hisat2_run_aurora";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}

my $shfile;
foreach ( 'AbsolutelyEmpty.sra', 'AbsolutelyEmpty.fastq',
	'AbsolutelyEmpty2.fastq' )
{
	push( @files, "$plugin_path/data/$_" );
	ok( -f "$plugin_path/data/$_", "infile $_ exists" );
	$shfile = "$plugin_path/data/$_";
	$shfile =~ s/\..*$/.sh/;
	unlink($shfile) if ( -f $shfile );
}
$coverage = "$plugin_path/data/AbsolutelyEmpty.chrlength";
ok( -f $coverage, "coverage input file" );
$genome = "$plugin_path/data";

$bigwigTracks = "$outpath/BigWig.html";
ok( !-f $bigwigTracks, "outfile does not exists before cmd run" );

@options = ( "A", "NOprojAtall" , 'project','dell');

$max_jobs = 40;

my $cmd =
    "perl -I $plugin_path/../lib  $exec "
  . " -max_jobs "
  . $max_jobs
 # . " -paired "    # or not?
  . " -files "
  . join( ' ', $files[0] )
  . " -bigwigTracks "
  . $bigwigTracks
  . " -sra "       # or not?
  . " -options "
  . join( ' ', @options )
  . " -coverage "
  . $coverage
  . " -genome "
  . $genome
  . " -outpath "
  . $outpath
  . " -debug";     ## never take this away!

system($cmd );

## this should have produced a script file

ok( -f "$plugin_path/data/AbsolutelyEmpty.sh", "script file created" );

@values =
  &read_file( "$plugin_path/data/AbsolutelyEmpty.sh", "script file created" );

#print "\$exp = " . root->print_perl_var_def( \@values ) . ";\n";
$exp = [
	'#! /bin/bash',
	'#SBATCH -n 10',
	'#SBATCH -N 1',
	'#SBATCH -t 02:00:00',
	'#SBATCH -A NOprojAtall',
	'#SBATCH -J AbsolutelyEmpty',
	'#SBATCH -o AbsolutelyEmpty%j.out',
	'#SBATCH -e AbsolutelyEmpty%j.err',
	'module load icc/2016.1.150-GCC-4.9.3-2.25 impi/5.1.2.150 SAMtools/1.3.1 '
		.'HISAT2/2.0.4 BEDTools/2.25.0 ucsc-tools/R2016a ',
	"hisat2  -x $plugin_path/data --sra-acc $plugin_path/data/AbsolutelyEmpty.sra --threads 10 "
		."--add-chrname > \$SNIC_TMP/AbsolutelyEmpty_hisat.sam",
	"samtools view -Sb  \$SNIC_TMP/AbsolutelyEmpty_hisat.sam | samtools sort -@ 9 -o \$SNIC_TMP/AbsolutelyEmpty_hisat.sorted.bam -",
	"if  [ -f \$SNIC_TMP/AbsolutelyEmpty_hisat.sorted.bam ]&&[ -s \$SNIC_TMP/AbsolutelyEmpty_hisat.sorted.bam ]; then",
	"rm -f \$SNIC_TMP/AbsolutelyEmpty_hisat.sam",
	'fi',
	"mv \$SNIC_TMP/AbsolutelyEmpty_hisat.sorted.bam $outpath/AbsolutelyEmpty_hisat.sorted.bam",
	"module purge",
	'module load GCC/4.9.3-2.25 OpenMPI/1.10.2 icc/2016.1.150-GCC-4.9.3-2.25 impi/5.1.2.150 BEDTools/2.25.0 Java/1.8.0_72 picard/2.8.2 ucsc-tools/R2016a ',
	"bedtools genomecov -bg -split -ibam $outpath/AbsolutelyEmpty_hisat.sorted.bam -g $plugin_path/data/AbsolutelyEmpty.chrlength |"
  		." sort -k1,1 -k2,2n > $outpath/AbsolutelyEmpty_hisat.bedGraph",
  	"polish_bed_like_files.pl -bed_file $outpath/AbsolutelyEmpty_hisat.bedGraph",
  	"bedGraphToBigWig $outpath/AbsolutelyEmpty_hisat.bedGraph $plugin_path/data/AbsolutelyEmpty.chrlength $outpath/AbsolutelyEmpty_hisat.bw",
  	
''
];
is_deeply( \@values, $exp, "right script for mapping one sra file" );

#print "\$exp = ".root->print_perl_var_def($value ).";\n";



sub hostname {
	open ( H, "hostname |" ) or die "I could not run the hostnale command\n$!\n";
	my $ret = join("",<H> );
	close ( H );
	$ret =~ s/\n//g;
	return $ret;
}

sub read_file {
	my $file    = shift;
	my $message = shift;
	my @ret;
	ok( -f $file, $message );
	if ( -f $file ) {
		open( IN, "<$file" ) or die $!;
		@ret = map { chomp; $_ } <IN>;
		close(IN);
	}
	return @ret;
}

