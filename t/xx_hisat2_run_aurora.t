#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 2;
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

@options = ( "A", "NOprojAtall" );

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
	"hisat2 -x $plugin_path/data --sra-acc $plugin_path/data/AbsolutelyEmpty.sra --threads 10 --add-chrname > $outpath/AbsolutelyEmpty_hisat.sam",
"samtools view -Sb  $outpath/AbsolutelyEmpty_hisat.sam | samtools sort -@ 9 -o $outpath/AbsolutelyEmpty_hisat.sorted.bam -",
"if  [ -f $outpath/AbsolutelyEmpty_hisat.sorted.bam ]&&[ -s $outpath/AbsolutelyEmpty_hisat.sorted.bam ]; then",
"rm -f $outpath/AbsolutelyEmpty_hisat.sam",
	'fi',
"bedtools genomecov -bg -split -ibam $outpath/AbsolutelyEmpty_hisat.sorted.bam -g $plugin_path/data/AbsolutelyEmpty.chrlength |"
  . " sort -k1,1 -k2,2n > $outpath/AbsolutelyEmpty_hisat.bedGraph",
	''
];
is_deeply( \@values, $exp, "right script for mapping one sra file" );

#print "\$exp = ".root->print_perl_var_def($value ).";\n";

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

