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
  . "-options n 1 N 1 t '00:02:00' p 2 proc 2 A 'snic2016-4-13' "
  . "-genome $plugin_path/data/hg38/hg38 "    ## does not even exist
  . "-coverage  $plugin_path/data/fake_hg38.chrom.sizes.txt "
  . "-bigwigTracks $plugin_path/data/output/hisat_BigWigTrackInfo.txt -dropDuplicates  -debug"
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
'track type=bigWig name="test_empty.fastq_hisat.sorted_picard_deduplicated" description="test_empty.fastq_hisat.sorted_picard_deduplicated" bigDataUrl=http://bone.bmc.lu.se/Public/test_empty.fastq_hisat.sorted_picard_deduplicated.bw'
];
is_deeply( $value, $exp, "The bigwig outfile contains the right entries" );

open( IN, "<$plugin_path/data/test_empty.fastq.sh" )
  or die "I could not read the created script file\n$!\n";
$value = [ map { chomp; $_ } <IN> ];
close(IN);
#print " \$exp = " . root->print_perl_var_def($value) . ";\n ";
$exp = [
	'#! /bin/bash',
	'#SBATCH -n 1',
	'#SBATCH -N 1',
	'#SBATCH -t 00:02:00',
	'#SBATCH -A snic2016-4-13',
	'#SBATCH -p 2',
	'#SBATCH -J test_empty.fastq',
	'#SBATCH -o test_empty.fastq.%j.out',
	'#SBATCH -e test_empty.fastq.%j.err',
	'module purge',
	'module load icc/2017.1.132-GCC-6.3.0-2.27 impi/2017.1.132 SAMtools/1.4.1 HISAT2/2.1.0 BEDTools/2.26.0 Java/1.8.0_92 ',
"hisat2  -x $plugin_path/data/hg38/hg38 -U $plugin_path/data/test_empty.fastq.gz --threads 2 --add-chrname > \$SNIC_TMP/test_empty.fastq_hisat.sam",
"samtools view -Sb  \$SNIC_TMP/test_empty.fastq_hisat.sam | samtools sort -@ 1 -o \$SNIC_TMP/test_empty.fastq_hisat.sorted.bam -",
"if  [ -f \$SNIC_TMP/test_empty.fastq_hisat.sorted.bam ]&&[ -s \$SNIC_TMP/test_empty.fastq_hisat.sorted.bam ]; then",
"rm -f \$SNIC_TMP/test_empty.fastq_hisat.sam",
	'fi',
	'module purge',
	'module load GCC/4.9.3-2.25 OpenMPI/1.10.2 icc/2016.1.150-GCC-4.9.3-2.25 impi/5.1.2.150 BEDTools/2.25.0 Java/1.8.0_72 picard/2.8.2 ucsc-tools/R2016a ',
	"picard MarkDuplicates INPUT=\$SNIC_TMP/test_empty.fastq_hisat.sorted.bam OUTPUT=\$SNIC_TMP/test_empty.fastq_hisat.sorted_picard_deduplicated.bam REMOVE_DUPLICATES=TRUE M=\$SNIC_TMP/test_empty.fastq_hisat.sorted_picard_metrix.txt",
"bedtools genomecov -bg -split -ibam \$SNIC_TMP/test_empty.fastq_hisat.sorted_picard_deduplicated.bam -g $plugin_path/data/fake_hg38.chrom.sizes.txt | sort -k1,1 -k2,2n > \$SNIC_TMP/test_empty.fastq_hisat.sorted_picard_deduplicated.bedGraph",
"polish_bed_like_files.pl -bed_file \$SNIC_TMP/test_empty.fastq_hisat.sorted_picard_deduplicated.bedGraph",
"bedGraphToBigWig \$SNIC_TMP/test_empty.fastq_hisat.sorted_picard_deduplicated.bedGraph /projects/fs1/med-sal/git/SLURM_bioinformatics_scripts/t/data/fake_hg38.chrom.sizes.txt \$SNIC_TMP/test_empty.fastq_hisat.sorted_picard_deduplicated.bw",
#"bedGraphToBigWig $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.bedGraph $plugin_path/data/fake_hg38.chrom.sizes.txt $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.bw",
"mv \$SNIC_TMP/test_empty.fastq_hisat.sorted.bam $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.sorted.bam",
"mv \$SNIC_TMP/test_empty.fastq_hisat.sorted_picard_deduplicated.bam $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.sorted_picard_deduplicated.bam",
"mv \$SNIC_TMP/test_empty.fastq_hisat.sorted_picard_deduplicated.bw $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.sorted_picard_deduplicated.bw",
	''
  ];

is_deeply( $value, $exp, "The script file contains the right entries" );


my $tmp =
    "perl -I $plugin_path/../lib $exec "
  . "-files $plugin_path/data/test_empty.fastq.gz "
  . "-options n 1 N 1 t '00:02:00' p 2 proc 2 A 'snic2016-4-13' "
  . "-genome $plugin_path/data/hg38/hg38 "    ## does not even exist
  . "-coverage  $plugin_path/data/fake_hg38.chrom.sizes.txt "
  . "-bigwigTracks $plugin_path/data/output/hisat_BigWigTrackInfo.txt -justMapping  -debug"
  . " > /dev/null"
  ;

system($tmp );

open( IN, "<$plugin_path/data/test_empty.fastq.sh" )
  or die "I could not read the created script file\n$!\n";
$value = [ map { chomp; $_ } <IN> ];
close(IN);

$tmp = [ @$exp[0..17,22], '' ];


is_deeply( $value, $tmp, "The script file contains the right entries -justMapping" );

#$exp = $value; ## I want to check the comments

#system ( "touch $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.sam");
#
#&run_cmd_and_read_script();
#@$exp[11] = "#@$exp[11]";
#
#is_deeply( $value, $exp, "The sam creation was blocked" );
#
#unlink( "$plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.sam" );
system ( "touch $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.sorted.bam");

&run_cmd_and_read_script(); ## all sam and sorted.bam creation should be blocked
foreach ( 11..15, 22) {
	@$exp[$_] = "#".@$exp[$_];
}
$exp = [ @$exp[0..15],
"cp $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.sorted.bam \$SNIC_TMP/test_empty.fastq_hisat.sorted.bam","",
 @$exp[16..(scalar(@$exp)-1)]
];

is_deeply( $value, $exp, "All sam and sorted.bam creation was blocked" );

#unlink( "$plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.sorted.bam");
system( "touch $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.sorted_picard_deduplicated.bam");
&run_cmd_and_read_script(); ## all sam and sorted.bam creation should be blocked
$exp = [ @$exp[0..20],
"cp $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.sorted_picard_deduplicated.bam \$SNIC_TMP/test_empty.fastq_hisat.sorted_picard_deduplicated.bam","",
 @$exp[21..(scalar(@$exp)-1)]
];

foreach (  20 , 27) {
	@$exp[$_] = "#".@$exp[$_];
}


is_deeply( $value, $exp, "All sam, sorted.bam and bedGraph creation was blocked" );

system( "touch $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.sorted_picard_deduplicated.bedGraph");
system( "touch $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.bw");

&run_cmd_and_read_script();

$exp = [ @$exp[0..24],
"cp $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.sorted_picard_deduplicated.bedGraph \$SNIC_TMP/test_empty.fastq_hisat.sorted_picard_deduplicated.bedGraph","",
 @$exp[25..(scalar(@$exp)-1)]
];

foreach (  23, 30) {
	@$exp[$_] = "#".@$exp[$_];
}

is_deeply( $value, $exp, "All sam, sorted.bam, bedGraph and bw creation was blocked" );

## and now a final check that the script is not run:

system( "touch $plugin_path/data/HISAT2_OUT/test_empty.fastq_hisat.sorted_picard_deduplicated.bw");
#$cmd =~ s/\-debug//;
&run_cmd_and_read_script(); ## all sam and sorted.bam creation should be blocked
for ( my $i = 0; $i <  @$exp;$i ++) {
	unless ( @$exp[$i] =~ m/^#/ ) {
		@$exp[$i] = "#@$exp[$i]";
	}
}
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