#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 4;

use FindBin;
my $plugin_path = $FindBin::Bin;
my ( $value, $exp );
my $exec = $plugin_path . "/../bin/bowtie2_run_aurora.pl";
ok( -f $exec, "the script has been found" );

if ( -d "$plugin_path/data/bowtie2/" ) {
	system( 'rm -Rf ' . "$plugin_path/data/bowtie2/" );
}
if ( -f -f "$plugin_path/data/test_empty.fastq.sh" ) {
	unlink( -f "$plugin_path/data/test_empty.fastq.sh" );
}
$value =
    "perl -I $plugin_path/../lib $exec "
  . "-files $plugin_path/data/test_empty.fastq.gz "
  . "-options n 1 N 1 t '00:02:00' p 2 "
  . "-genome $plugin_path/data/hg38/hg38 "    ## does not even exist
  . "-coverage  $plugin_path/data/fake_hg38.chrom.sizes.txt "
  . "-bigwigTracks $plugin_path/data/output/BigWigTrackInfo.txt  -debug";
print "the command:\n$value\n";
system($value );

ok( -d "$plugin_path/data/bowtie2/", "the script created the outpath" );

ok( -f "$plugin_path/data/test_empty.fastq.sh",
	"the SBATCH script was created" );

open( IN, "<$plugin_path/data/test_empty.fastq.sh" );
$value = [ map { chomp; $_ } <IN> ];
close(IN);
$exp = [
	'#! /bin/bash',
	'#SBATCH -n 1',
	'#SBATCH -N 1',
	'#SBATCH -t 00:02:00',
	'#SBATCH -J test_empty.fastq',
	'#SBATCH -o test_empty.fastq%j.out',
	'#SBATCH -e test_empty.fastq%j.err',
"bowtie2 -x $plugin_path/data/hg38/hg38  -p 2 -U '$plugin_path/data/test_empty.fastq.gz' -S $plugin_path/data/bowtie2/test_empty.fastq_bowtie2.sam",
"samtools view -Sb  $plugin_path/data/bowtie2/test_empty.fastq_bowtie2.sam | samtools sort -@ 1 - $plugin_path/data/bowtie2/test_empty.fastq_bowtie2.sorted",
"if  [ -f $plugin_path/data/bowtie2/test_empty.fastq_bowtie2.sorted.bam ]&&[ -s $plugin_path/data/bowtie2/test_empty.fastq_bowtie2.sorted.bam ]; then",
"rm -f $plugin_path/data/bowtie2/test_empty.fastq_bowtie2.sam",
	'fi',
	"bedtools genomecov -bg -split -ibam $plugin_path/data/bowtie2/test_empty.fastq_bowtie2.sorted.bam -g $plugin_path/data/fake_hg38.chrom.sizes.txt | sort -k1,1 -k2,2n > $plugin_path/data/bowtie2/test_empty.fastq_bowtie2.bedGraph",
	"bedGraphToBigWig $plugin_path/data/bowtie2/test_empty.fastq_bowtie2.bedGraph $plugin_path/data/fake_hg38.chrom.sizes.txt $plugin_path/data/bowtie2/test_empty.fastq_bowtie2.bw",
	''
];

#print " \$exp = " . root->print_perl_var_def($value) . "\n ";

is_deeply( $value, $exp, "The script contains the right entries" )



