#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 10;

use FindBin;
my $plugin_path = $FindBin::Bin;
my ( $value, $exp, $cmd );
my $exec = $plugin_path . "/../bin/bowtie_run_aurora.pl";
ok( -f $exec, "the script has been found" );

if ( -d "$plugin_path/data/bowtie/" ) {
	system( 'rm -Rf ' . "$plugin_path/data/bowtie/" );
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
  . "-bigwigTracks $plugin_path/data/output/BigWigTrackInfo.txt  -debug > /dev/null 2> /dev/null";
print "the command:\n$cmd\n";
system($cmd );

ok( -d "$plugin_path/data/bowtie/", "the script created the outpath" );

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
	"module purge",
	 'module load icc/2015.3.187-GNU-4.9.3-2.25 impi/5.0.3.048 Bowtie/1.1.2 SAMtools/0.1.20 ',
"bowtie  -p 2 --sam $plugin_path/data/hg38/hg38 '$plugin_path/data/test_empty.fastq.gz' > $plugin_path/data/bowtie/test_empty.fastq_bowtie.sam",
"samtools view -Sb  $plugin_path/data/bowtie/test_empty.fastq_bowtie.sam | samtools sort -@ 1 -o $plugin_path/data/bowtie/test_empty.fastq_bowtie.sorted.bam -",
"if  [ -f $plugin_path/data/bowtie/test_empty.fastq_bowtie.sorted.bam ]&&[ -s $plugin_path/data/bowtie/test_empty.fastq_bowtie.sorted.bam ]; then",
"rm -f $plugin_path/data/bowtie/test_empty.fastq_bowtie.sam",
	'fi','',
	"#bedtools genomecov -bg -split -ibam $plugin_path/data/bowtie/test_empty.fastq_bowtie.sorted.bam -g $plugin_path/data/fake_hg38.chrom.sizes.txt | sort -k1,1 -k2,2n > $plugin_path/data/bowtie/test_empty.fastq_bowtie.bedGraph",
	"#polish_bed_like_files.pl -bed_file $plugin_path/data/bowtie/test_empty.fastq_bowtie.bedGraph",
	"#bedGraphToBigWig $plugin_path/data/bowtie/test_empty.fastq_bowtie.bedGraph $plugin_path/data/fake_hg38.chrom.sizes.txt $plugin_path/data/bowtie/test_empty.fastq_bowtie.bw"
];

is_deeply( $value, $exp, "The script contains the right entries" );


open ( IN, "<$plugin_path/data/output/BigWigTrackInfo.txt" ) or die "I could not read the created BigWigTrackInfo file\n$!\n";
$value = [ map { chomp; $_ } <IN> ];
close(IN);

#die " \$exp = " . root->print_perl_var_def($value) . "\n ";

my $exp2 = ['track type=bigWig name="test_empty.fastq_bowtie" description="test_empty.fastq_bowtie" bigDataUrl=http://bone.bmc.lu.se/Public/test_empty.fastq_bowtie.bw' ];
is_deeply( $value, $exp2, "The bigwig outfile contains the right entries" );

#### now lets check the script if a outfile is already present 

system( "touch $plugin_path/data/bowtie/test_empty.fastq_bowtie.sam" );
system($cmd );

open( IN, "<$plugin_path/data/test_empty.fastq.sh" );
$value = [ map { chomp; $_ } <IN> ];
close(IN);

@$exp[10] = "#@$exp[10]";

is_deeply( $value, $exp, "The script contains the right entries ( .sam existing)" );


system( "touch $plugin_path/data/bowtie/test_empty.fastq_bowtie.sorted.bam" );
unlink( "$plugin_path/data/bowtie/test_empty.fastq_bowtie.sam" );
system($cmd );

open( IN, "<$plugin_path/data/test_empty.fastq.sh" );
$value = [ map { chomp; $_ } <IN> ];
close(IN);
foreach ( 11..14) {
	@$exp[$_] = "#@$exp[$_]";
}
is_deeply( $value, $exp, "The script contains the right entries ( .sorted.bam existing)" );


system( "touch $plugin_path/data/bowtie/test_empty.fastq_bowtie.bedGraph" );
system($cmd );

open( IN, "<$plugin_path/data/test_empty.fastq.sh" );
$value = [ map { chomp; $_ } <IN> ];
close(IN);
@$exp[16] = "#@$exp[16]";


is_deeply( $value, $exp, "The script contains the right entries ( .bedGraph existing)" );


system( "touch $plugin_path/data/bowtie/test_empty.fastq_bowtie.bw" );
system($cmd );

open( IN, "<$plugin_path/data/test_empty.fastq.sh" );
$value = [ map { chomp; $_ } <IN> ];
close(IN);
foreach ( 17,18 ) {
	@$exp[$_] = "#@$exp[$_]";
}
is_deeply( $value, $exp, "The script contains the right entries ( .bw existing)" );


$cmd =
    "perl -I $plugin_path/../lib $exec "
  . "-files $plugin_path/data/test_empty.fastq.gz "
  . "-options n 1 N 1 t '00:02:00' p 2 A 'snic2016-4-13' "
  . "-bowtie_options '-m 500 -v2 -best -strata' "
  . "-genome $plugin_path/data/hg38/hg38 "    ## does not even exist
  . "-coverage  $plugin_path/data/fake_hg38.chrom.sizes.txt "
  . "-bigwigTracks $plugin_path/data/output/BigWigTrackInfo.txt  -debug > /dev/null 2>/dev/null";
print "the command:\n$cmd\n";
system($cmd );


open( IN, "<$plugin_path/data/test_empty.fastq.sh" );
$value = [ map { chomp; $_ } <IN> ];
close(IN);

@$exp[10] = "#bowtie -m 500 -v2 -best -strata -p 2 --sam $plugin_path/data/hg38/hg38 '$plugin_path/data/test_empty.fastq.gz' > $plugin_path/data/bowtie/test_empty.fastq_bowtie.sam",

is_deeply( $value, $exp, "extra bowtie_options end up in right place" );


#print " \$exp = " . root->print_perl_var_def($value) . "\n ";

