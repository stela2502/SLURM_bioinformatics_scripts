#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 11;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my (
	$value,        @values,  $exp,      $max_jobs, @files,
	$bigwigTracks, @options, $coverage, $genome,   $outpath,
);

my $exec = $plugin_path . "/../bin/star_run_aurora.pl";
ok( -f $exec, 'the script has been found' );
$outpath = "$plugin_path/data/output/star_run_aurora";
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
  . " -justMapping -max_jobs "
  . $max_jobs
 # . " -paired "    # or not?
  . " -files "
  . join( ' ', $files[0] )
  . " -bigwigTracks "
  . $bigwigTracks
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
 '#SBATCH -J Ab42ce', 
 '#SBATCH -o Ab42ce.%j.out', 
 '#SBATCH -e Ab42ce.%j.err', 
 'module purge', 
 'module load GCC/5.4.0-2.26 OpenMPI/1.10.3 STAR/2.6.0c SAMtools/1.4 BEDTools/2.26.0 ', 
 'date', 
 "STAR  --genomeDir $plugin_path/data --readFilesIn $plugin_path/data/AbsolutelyEmpty.sra --runThreadN 10  --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMmultNmax 1 --outFileNamePrefix \$SNIC_TMP/AbsolutelyEmpty --outSAMtype BAM SortedByCoordinate --twopassMode Basic ", 
 'date', 'date', 
 '## here comes a subscript entry:', 
 'module purge', 
 'module load GCC/4.9.3-2.25 OpenMPI/1.10.2 icc/2016.1.150-GCC-4.9.3-2.25 impi/5.1.2.150 BEDTools/2.25.0 Java/1.8.0_72 picard/2.8.2 ucsc-tools/R2016a ', 
 'if [ -f $SNIC_TMP/\'AbsolutelyEmptyAligned.sortedByCoord.out.bam\' ];then ', 
 "  cp \$SNIC_TMP/AbsolutelyEmptyAligned.sortedByCoord.out.bam $plugin_path/data/output/star_run_aurora/AbsolutelyEmptyAligned.sortedByCoord.out.bam", 
 'else', 
 '  ( >&2 echo \'file $SNIC_TMP/AbsolutelyEmptyAligned.sortedByCoord.out.bam was not produced - files in the same path:\')', '  for filename in $SNIC_TMP/*',
 '  do', 
 '    (>&2 echo $filename) ', 
 '  done;', 
 'fi', 
 'date', 
 'exit 0'
];
is_deeply( \@values, $exp, "right script for mapping one sra file" );

#print "\$exp = ".root->print_perl_var_def($value ).";\n";

## now lets create the main outfile and check the script again...
open ( OUT , ">$plugin_path/data/output/star_run_aurora/AbsolutelyEmptyAligned.sortedByCoord.out.bam") or die $!;
print OUT "Just a test";
close ( OUT );

system( $cmd );

@values =
  &read_file( "$plugin_path/data/AbsolutelyEmpty.sh", "script file created #2" );

#print "\$exp = " . root->print_perl_var_def( \@values ) . ";\n";

$exp = [
 '#! /bin/bash', 
 '#SBATCH -n 10', 
 '#SBATCH -N 1', 
 '#SBATCH -t 02:00:00', 
 '#SBATCH -A NOprojAtall', 
 '#SBATCH -J Ab42ce', 
 '#SBATCH -o Ab42ce.%j.out', 
 '#SBATCH -e Ab42ce.%j.err', 
 'module purge', 
 'module load GCC/5.4.0-2.26 OpenMPI/1.10.3 STAR/2.6.0c SAMtools/1.4 BEDTools/2.26.0 ', 
 '#date', 
 "#STAR  --genomeDir $plugin_path/data --readFilesIn $plugin_path/data/AbsolutelyEmpty.sra --runThreadN 10  --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMmultNmax 1 --outFileNamePrefix \$SNIC_TMP/AbsolutelyEmpty --outSAMtype BAM SortedByCoordinate --twopassMode Basic ", 
 '#date', 'date', 
 '## here comes a subscript entry:', 
 'module purge', 
 'module load GCC/4.9.3-2.25 OpenMPI/1.10.2 icc/2016.1.150-GCC-4.9.3-2.25 impi/5.1.2.150 BEDTools/2.25.0 Java/1.8.0_72 picard/2.8.2 ucsc-tools/R2016a ', 
 '#if [ -f $SNIC_TMP/\'AbsolutelyEmptyAligned.sortedByCoord.out.bam\' ];then ', 
 "#  cp \$SNIC_TMP/AbsolutelyEmptyAligned.sortedByCoord.out.bam $plugin_path/data/output/star_run_aurora/AbsolutelyEmptyAligned.sortedByCoord.out.bam", 
 '#else', 
 '#  ( >&2 echo \'file $SNIC_TMP/AbsolutelyEmptyAligned.sortedByCoord.out.bam was not produced - files in the same path:\')', 
 '#  for filename in $SNIC_TMP/*',
 '#  do', 
 '#    (>&2 echo $filename) ', 
 '#  done;', 
 '#fi', 
 'date', 
 'exit 0'
];
is_deeply( \@values, $exp, "Do nothing script if outfile exists" );

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

