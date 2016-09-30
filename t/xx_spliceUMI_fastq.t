#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 42;
use stefans_libs::flexible_data_structures::data_table;
use Digest::MD5;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $data_read, $sampleDefinition, $UMI_read, @options, );

my $exec = $plugin_path . "/../bin/spliceUMI_fastq.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/spliceUMI_fastq";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}
$data_read = $plugin_path."/data/spliceUMI/GMLP_S1_R2_001.part.fastq";
$UMI_read = $plugin_path."/data/spliceUMI/GMLP_S1_R1_001.part.fastq";
$sampleDefinition = $plugin_path."/data/spliceUMI/sample_sheet_160516.csv";

ok ( -f $data_read, "data reads file '$data_read'");
ok ( -f $UMI_read, "UMI reads file");
ok ( -f $sampleDefinition, "sample definition file");

my $cmd =
    "perl -I $plugin_path/../lib  $exec "
. " -outpath " . $outpath 
. " -data_read " . $data_read 
. " -sampleDefinition " . $sampleDefinition 
. " -UMI_read " . $UMI_read 
. " -options " . join(' ', @options )
. " -debug";

system ( $cmd ) ;
my @digests = (
'18ebb51714e774d507cf986d22973b59',  'GMLP_S1_R2_001.part.AGATAG.fastq',
'84e22846c897b926a076e3f72e6b9750',  'GMLP_S1_R2_001.part.ATACGA.fastq',
'18b784a28bac7813e5ac2c63a7f0a766',  'GMLP_S1_R2_001.part.CAGATA.fastq',
'09e7742d15f941119f73a80d40e23f91',  'GMLP_S1_R2_001.part.CCTGAG.fastq',
'1595a98b81b3abdc521b791a4e448f24',  'GMLP_S1_R2_001.part.GAGTCC.fastq',
'a430870b56ae8ade84aaa632f2aa58f8',  'GMLP_S1_R2_001.part.GTCGCT.fastq',
'7728f8a63c3bb602263decd7e921788c',  'GMLP_S1_R2_001.part.NoSample.fastq',
'16a5f97a26f86c01ca9cc029b714551b',  'GMLP_S1_R2_001.part.TCTCGC.fastq',
'd4f56cccaaa9f8757f88925dbdd53c47',  'GMLP_S1_R2_001.part.TGCATT.fastq'
);

foreach ( qw(  GMLP_S1_R2_001.part.CCTGAG.fastq    GMLP_S1_R2_001.part.TCTCGC.fastq
GMLP_S1_R2_001.part.AGATAG.fastq  GMLP_S1_R2_001.part.GAGTCC.fastq    GMLP_S1_R2_001.part.TGCATT.fastq
GMLP_S1_R2_001.part.ATACGA.fastq  GMLP_S1_R2_001.part.GTCGCT.fastq
GMLP_S1_R2_001.part.CAGATA.fastq  GMLP_S1_R2_001.part.NoSample.fastq ) ) {
	ok ( -f $outpath."/".$_, "outfile $_" );
}

my ($ctx, $filehandle, $tmp);
for ( my $i = 0; $i < @digests; $i +=2 ) {
	$ctx = Digest::MD5->new;
	open ( $filehandle, "<$outpath/$digests[$i+1]" );
	$ctx->addfile($filehandle);
	close ( $filehandle );
	ok( $digests[$i] eq  ($tmp = $ctx->hexdigest), "outfile $digests[$i+1] ($tmp)");
}


$data_read .= ".gz";
$UMI_read  .= ".gz";
ok ( -f $data_read, "gzipped data reads file '$data_read'");
ok ( -f $UMI_read, "gzipped UMI reads file");

$cmd =
    "perl -I $plugin_path/../lib  $exec "
. " -outpath " . $outpath 
. " -data_read " . $data_read 
. " -sampleDefinition " . $sampleDefinition 
. " -UMI_read " . $UMI_read 
. " -options " . join(' ', @options )
. " -debug";

system ( $cmd ) ;

foreach ( qw(  GMLP_S1_R2_001.part.CCTGAG.fastq    GMLP_S1_R2_001.part.TCTCGC.fastq
GMLP_S1_R2_001.part.AGATAG.fastq  GMLP_S1_R2_001.part.GAGTCC.fastq    GMLP_S1_R2_001.part.TGCATT.fastq
GMLP_S1_R2_001.part.ATACGA.fastq  GMLP_S1_R2_001.part.GTCGCT.fastq
GMLP_S1_R2_001.part.CAGATA.fastq  GMLP_S1_R2_001.part.NoSample.fastq ) ) {
	ok ( -f $outpath."/".$_.".gz", "outfile $_.gz" );
}


for ( my $i = 0; $i < @digests; $i +=2 ) {
	$ctx = Digest::MD5->new;
	open ( $filehandle, "zcat $outpath/$digests[$i+1].gz |" );
	$ctx->addfile($filehandle);
	close ( $filehandle );
	ok( $digests[$i] eq  ($tmp = $ctx->hexdigest), "unzipped gzipped outfile $digests[$i+1].gz ($tmp)");
}

#print "\$exp = ".root->print_perl_var_def($value ).";\n