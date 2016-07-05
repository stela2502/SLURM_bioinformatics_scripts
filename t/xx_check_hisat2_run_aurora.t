#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 9;

use FindBin;
my $plugin_path = $FindBin::Bin;
my ( $value, $exp, $cmd );
my $exec = $plugin_path . "/../bin/hisat2_run_aurora.pl";
ok( -f $exec, "the script has been found" );

if ( -d "$plugin_path/data/HISAT2_OUT/" ) {
	system( 'rm -Rf ' . "$plugin_path/data/HISAT2_OUT/" );
}
if ( -f "$plugin_path/data/test_empty.fastq.sh" ) {
	unlink( -f "$plugin_path/data/test_empty.fastq.sh" );
}
$cmd =
    "perl -I $plugin_path/../lib $exec "
  . "-files $plugin_path/data/test_empty.fastq.gz "
  . "-options n 1 N 1 t '00:02:00' p 2 proc 2 "
  . "-genome $plugin_path/data/hg38/hg38 "    ## does not even exist
  . "-coverage  $plugin_path/data/fake_hg38.chrom.sizes.txt "
  . "-bigwigTracks $plugin_path/data/output/hisat_BigWigTrackInfo.txt  -debug"
  . " > /dev/null"
  ;
print "the command:\n$cmd\n";

system($cmd );

ok( -d "$plugin_path/data/HISAT2_OUT/", "outpath created" );

open( IN, "<$plugin_path/data/output/hisat_BigWigTrackInfo.txt" )
  or die "I could not read the created BigWigTrackInfo file\n$!\n";
$value = [ map { chomp; $_ } <IN> ];
close(IN);

#print " \$exp = " . root->print_perl_var_def($value) . ";\n ";
$exp = [
'track type=bigWig name="test_empty.fastq_hisat" description="test_empty.fastq_hisat" bigDataUrl=http://bone.bmc.lu.se/Public/test_empty.fastq_hisat.bw'
];
is_deeply( $value, $exp, "The bigwig outfile contains the right entries" );

open( IN, "<$plugin_path/data/test_empty.fastq.sh" )
  or die "I could not read the created script file\n$!\n";
$value = [ map { chomp; $_ } <IN> ];
close(IN);
print " \$exp = " . root->print_perl_var_def($value) . ";\n ";
$exp = [
	'#! /bin/bash',
	'#SBATCH -n 1',
	'#SBATCH -N 1',
	'#SBATCH -t 00:02:00',
	'#SBATCH -J test_empty.fastq',
	'#SBATCH -o test_empty.fastq%j.out',
	'#SBATCH -e test_empty.fastq%j.err',
"hisat2 -x $plugin_path/data/hg38/hg38 -U $plugin_path/data/test_empty.fastq.gz --threads 2 --add-chrname > $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.sam",
"samtools view -Sb  $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.sam | samtools sort -@ 1 -o $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.sorted.bam -",
"if  [ -f $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.sorted.bam ]&&[ -s $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.sorted.bam ]; then",
"rm -f $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.sam",
	'fi',
"bedtools genomecov -bg -split -ibam $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.sorted.bam -g $plugin_path/data/fake_hg38.chrom.sizes.txt | sort -k1,1 -k2,2n > $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.bedGraph",
'polish_bed_like_files.pl -bed_file /home/med-sal/git_Projects/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.bedGraph',
"bedGraphToBigWig $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.bedGraph $plugin_path/data/fake_hg38.chrom.sizes.txt $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.bw",
	''
  ];

is_deeply( $value, $exp, "The script file contains the right entries" );

$exp = $value; ## I want to check the comments

system ( "touch $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.sam");

&run_cmd_and_read_script();
@$exp[7] = "#@$exp[7]";

is_deeply( $value, $exp, "The sam creation was blocked" );

unlink( "$plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.sam" );
system ( "touch $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.sorted.bam");

&run_cmd_and_read_script(); ## all sam and sorted.bam creation should be blocked
@$exp[8] = "#@$exp[8]";
@$exp[9] = "#@$exp[9]";
@$exp[10] = "#@$exp[10]";
@$exp[11] = "#@$exp[11]";

is_deeply( $value, $exp, "All sam and sorted.bam creation was blocked" );

system( "touch $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.bedGraph");
&run_cmd_and_read_script(); ## all sam and sorted.bam creation should be blocked
@$exp[12] = "#@$exp[12]";

is_deeply( $value, $exp, "All sam, sorted.bam and bedGraph creation was blocked" );

system( "touch $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.bw");

&run_cmd_and_read_script();
@$exp[13] = "#@$exp[13]";
if ( defined @$exp[14] ) {
	@$exp[14] = "#@$exp[14]";
}
is_deeply( $value, $exp, "All sam, sorted.bam, bedGraph and bw creation was blocked" );

## and now a final check that the script is not run:

$cmd =~ s/\-debug//;
&run_cmd_and_read_script(); ## all sam and sorted.bam creation should be blocked
is_deeply( $value, $exp, "All sam, sorted.bam and bedGraph creation was blocked and you have not seen an SBATCH error message - right?" );

#unlink ( "$plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.bedGraph");
#&run_cmd_and_read_script(); ## all sam and sorted.bam creation should be blocked


sub run_cmd_and_read_script {
	system($cmd );

	open( IN, "<$plugin_path/data/test_empty.fastq.sh" )
  	or die "I could not read the created script file\n$!\n";
	$value = [ map { chomp; $_ } <IN> ];
	close(IN);
}