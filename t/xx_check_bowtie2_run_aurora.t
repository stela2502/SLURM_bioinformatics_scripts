#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 9;

use FindBin;
my $plugin_path = $FindBin::Bin;
my ( $value, $exp, $cmd );
my $exec = $plugin_path . "/../bin/bowtie2_run_aurora.pl";
ok( -f $exec, "the script has been found" );

if ( -d "$plugin_path/data/bowtie2/" ) {
	system( 'rm -Rf ' . "$plugin_path/data/bowtie2/" );
}
if ( -f  "$plugin_path/data/test_empty.fastq.sh" ) {
	unlink( -f "$plugin_path/data/test_empty.fastq.sh" );
}
$cmd =
    "perl -I $plugin_path/../lib $exec "
  . "-files $plugin_path/data/test_empty.fastq.gz "
  . "-options n 1 N 1 t '00:02:00' p 2 A 'snic2016-4-13' "
  . "-genome $plugin_path/data/hg38/hg38 "    ## does not even exist
  . "-coverage  $plugin_path/data/fake_hg38.chrom.sizes.txt "
  . "-bigwigTracks $plugin_path/data/output/BigWigTrackInfo.txt  -debug > /dev/null";
print "the command:\n$cmd\n";
system($cmd );

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
	'#SBATCH -A snic2016-4-13',
	'#SBATCH -J test_empty.fastq',
	'#SBATCH -o test_empty.fastq%j.out',
	'#SBATCH -e test_empty.fastq%j.err',
	'module load icc/2016.1.150-GCC-4.9.3-2.25 impi/5.1.2.150  Bowtie2/2.2.8 BEDTools/2.25.0 picard/2.8.2 ucsc-tools/R2016a '
	.'icc/2016.1.150-GCC-4.9.3-2.25 impi/5.1.2.150 HISAT2/2.0.4 BEDTools/2.25.0 ucsc-tools/R2016a ',
"bowtie2 -x $plugin_path/data/hg38/hg38  -p 2 -U '$plugin_path/data/test_empty.fastq.gz' -S $plugin_path/data/bowtie2/test_empty.fastq_bowtie2.sam",
"samtools view -Sb  $plugin_path/data/bowtie2/test_empty.fastq_bowtie2.sam | samtools sort -@ 1 -o $plugin_path/data/bowtie2/test_empty.fastq_bowtie2.sorted.bam -",
"if  [ -f $plugin_path/data/bowtie2/test_empty.fastq_bowtie2.sorted.bam ]&&[ -s $plugin_path/data/bowtie2/test_empty.fastq_bowtie2.sorted.bam ]; then",
"rm -f $plugin_path/data/bowtie2/test_empty.fastq_bowtie2.sam",
	'fi','',
	"bedtools genomecov -bg -split -ibam $plugin_path/data/bowtie2/test_empty.fastq_bowtie2.sorted.bam -g $plugin_path/data/fake_hg38.chrom.sizes.txt | sort -k1,1 -k2,2n > $plugin_path/data/bowtie2/test_empty.fastq_bowtie2.bedGraph",
	"polish_bed_like_files.pl -bed_file $plugin_path/data/bowtie2/test_empty.fastq_bowtie2.bedGraph",
	"bedGraphToBigWig $plugin_path/data/bowtie2/test_empty.fastq_bowtie2.bedGraph $plugin_path/data/fake_hg38.chrom.sizes.txt $plugin_path/data/bowtie2/test_empty.fastq_bowtie2.bw",
	''
];

is_deeply( $value, $exp, "The script contains the right entries" );


open ( IN, "<$plugin_path/data/output/BigWigTrackInfo.txt" ) or die "I could not read the created BigWigTrackInfo file\n$!\n";
$value = [ map { chomp; $_ } <IN> ];
close(IN);

#die " \$exp = " . root->print_perl_var_def($value) . "\n ";

my $exp2 = ['track type=bigWig name="test_empty.fastq_bowtie2" description="test_empty.fastq_bowtie2" bigDataUrl=http://bone.bmc.lu.se/Public/test_empty.fastq_bowtie2.bw' ];
is_deeply( $value, $exp2, "The bigwig outfile contains the right entries" );

#### now lets check the script if a outfile is already present 

system( "touch $plugin_path/data/bowtie2/test_empty.fastq_bowtie2.sam" );
system($cmd );

open( IN, "<$plugin_path/data/test_empty.fastq.sh" );
$value = [ map { chomp; $_ } <IN> ];
close(IN);

@$exp[9] = "#@$exp[9]";

is_deeply( $value, $exp, "The script contains the right entries ( .sam existing)" );


system( "touch $plugin_path/data/bowtie2/test_empty.fastq_bowtie2.sorted.bam" );
unlink( "$plugin_path/data/bowtie2/test_empty.fastq_bowtie2.sam" );
system($cmd );

open( IN, "<$plugin_path/data/test_empty.fastq.sh" );
$value = [ map { chomp; $_ } <IN> ];
close(IN);
foreach ( 10..13) {
	@$exp[$_] = "#@$exp[$_]";
}
is_deeply( $value, $exp, "The script contains the right entries ( .sorted.bam existing)" );


system( "touch $plugin_path/data/bowtie2/test_empty.fastq_bowtie2.bedGraph" );
system($cmd );

open( IN, "<$plugin_path/data/test_empty.fastq.sh" );
$value = [ map { chomp; $_ } <IN> ];
close(IN);
foreach ( 15 ) {
	@$exp[$_] = "#@$exp[$_]";
}
is_deeply( $value, $exp, "The script contains the right entries ( .bedGraph existing)" );


system( "touch $plugin_path/data/bowtie2/test_empty.fastq_bowtie2.bw" );
system($cmd );

open( IN, "<$plugin_path/data/test_empty.fastq.sh" );
$value = [ map { chomp; $_ } <IN> ];
close(IN);
foreach ( 16, 17 ) {
	@$exp[$_] = "#@$exp[$_]";
}
is_deeply( $value, $exp, "The script contains the right entries ( .bw existing)" );

#print " \$exp = " . root->print_perl_var_def($value) . "\n ";

