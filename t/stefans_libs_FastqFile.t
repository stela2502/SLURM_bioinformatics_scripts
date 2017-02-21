#! /usr/bin/perl
use strict;
use warnings;
use Test::More tests => 9;
BEGIN { use_ok 'stefans_libs::FastqFile' }
use stefans_libs::root;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp );
my $OBJ = stefans_libs::FastqFile->new( { 'debug' => 1 } );
is_deeply( ref($OBJ), 'stefans_libs::FastqFile',
	'simple test of function stefans_libs::FastqFile -> new() ' );

#print "\$exp = ".root->print_perl_var_def($value ).";\n";
$value = $OBJ->open_file( $plugin_path . "/data/repeat.fastq.gz" );
is_deeply( ref($value), 'GLOB', 'opened gzipped txt file' );
@values = <$value>;
close($value);

#print "\$exp = ".root->print_perl_var_def([@values[0..3]]).";\n";
$exp = [
	'@M04223:22:000000000-AUY26:1:1101:15371:1335 1:N:0:1
', 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
', '+
', '1>>>>>>1>>>>AAAAAA>>>>>>>><<<<B@@<<<<<<<:::99999999
'
];
is_deeply( [ @values[ 0 .. 3 ] ], $exp, "the right sequnce has been read\n" );

$value = $OBJ->open_file( $plugin_path . "/data/repeat.fastq" );

is_deeply( ref($value), 'GLOB', 'opened normal txt file' );
@values = <$value>;
close($value);
is_deeply( [ @values[ 0 .. 3 ] ], $exp, "the right sequnce has been read\n" );

$OBJ->exclude_read(
	$plugin_path . "/data/repeat.fastq.gz",
	'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT',
	$plugin_path . "/data/outpath/repeat.fastq"
);
$value  = $OBJ->open_file( $plugin_path . "/data/outpath/repeat.fastq" );
@values = <$value>;
close($value);

#print "\$exp = ".root->print_perl_var_def([@values]).";\n";
$exp = [
	'@M04223:22:000000000-AUY26:1:1101:16818:1622 1:N:0:1
', 'TCTTTTTTTTCTTTTTTTTTTTTTTTTTTTTTTTCTTCCTTTTTTCTTCTT
', '+
', '111>>11>>A013B11A0/A//A/A/>/>/>E//011211B211/0B12BF
', '@M04223:22:000000000-AUY26:1:1101:15584:1697 1:N:0:1
', 'TGTGGCACGATGTGACGAGGCGGCCGAGTGGTTAAGGCGATGGCGGTTCAG
', '+
', '>1>1111B11>1F13A1100A00AA///F//FFB111A//E//0BGAEGBB
', '@M04223:22:000000000-AUY26:1:1101:12954:1698 1:N:0:1
', 'TTCGGCAAGTCAAGCAGTTAACGCGGGTTCGATTCCCGGGAGATCGGAAGA
', '+
', '>A1111111BB@311A1FGEDE0A0?CEGF?AGHFDF???EECE////>EG
'
];

is_deeply( [@values], $exp, "the right sequnce has been read\n" );

## and now test the extract_barcode function

#TGTGGCACGATGTGACGAGGCGGCCGAGTGGTTAAGGCGATGGACTGCTAA

my $adapter        = 'AGATCGGAAGAGCGGTTCAG';
my $sample_barcode = 'GGCA';

$OBJ->extract_barcode( $plugin_path . "/data/repeat.fastq.gz",
	$adapter, $sample_barcode, 0,
	$plugin_path . "/data/outpath/Barcodes_GGCA.fastq" );
$value  = $OBJ->open_file( $plugin_path . "/data/outpath/Barcodes_GGCA.fastq" );
@values = <$value>;
close($value);

#print "\$exp = ".root->print_perl_var_def([@values]).";\n";
$exp = [
	'@M04223:22:000000000-AUY26:1:1101:15584:1697:TGTCG 1:N:0:1
', 'ATGTGACGAGGCGGCCGAGTGGTTAAGGCGATGGCGGTTCAG
', '+
', '1>1F13A1100A00AA///F//FFB111A//E//0BGAEGBB
', '@M04223:22:000000000-AUY26:1:1101:12954:1698:AD11:TTCAG 1:N:0:1
', 'TCAAGCAGTTAACGCGGGTTCGATTCCCGGG
', '+
', 'BB@311A1FGEDE0A0?CEGF?AGHFDF???
'
];

is_deeply( \@values, $exp, "right sequences returned" );

$OBJ->select_4_str( $plugin_path . "/data/outpath/Barcodes_GGCA.fastq",
	":AD\\d", 0, $plugin_path . "/data/outpath/Barcodes_GGCA_AD.fastq" );
$value =
  $OBJ->open_file( $plugin_path . "/data/outpath/Barcodes_GGCA_AD.fastq" );
@values = <$value>;
close($value);
$exp = [
	'@M04223:22:000000000-AUY26:1:1101:12954:1698:AD11:TTCAG 1:N:0:1
', 'TCAAGCAGTTAACGCGGGTTCGATTCCCGGG
', '+
', 'BB@311A1FGEDE0A0?CEGF?AGHFDF???
'
];
is_deeply( \@values, $exp, "right select_4_str sequences returned" );

my $result = [];
my $func   = sub {
	my ( $fastqfile, @entries ) = @_;
	for ( my $i = 0 ; $i < @entries ; $i++ ) {
		@$result[$i] ||= [];
		push( @{ @$result[$i] }, $entries[$i]->copy() );
		$entries[$i]->clear();
	}
};

ok( -f "$plugin_path/data/R1.fastq.gz", "$plugin_path/data/R1.fastq.gz" );
ok( -f "$plugin_path/data/R2.fastq.gz", "$plugin_path/data/R2.fastq.gz" );
ok( -f "$plugin_path/data/I1.fastq.gz", "$plugin_path/data/I1.fastq.gz" );

$OBJ->filter_multiple_files(
	$func,
	"$plugin_path/data/R1.fastq.gz",
	"$plugin_path/data/R2.fastq.gz",
	"$plugin_path/data/I1.fastq.gz"
);

ok( @{ @$result[0] } == 5, "right number of fastq entries" );

#print join( " ", @$result[0] ) . "\n";
$exp = [
	[
		'@NB501227:57:HGJM5BGX2:1:11101:3408:1053 1:N:0:ATGAATCT',
		'TGAGAGGGTGCTAGCCACGGAGTAAN', '+', 'AAAAAEEEEEAEEEEEEEEEEEEEE#'
	],
	[
		'@NB501227:57:HGJM5BGX2:1:11101:4228:1060 1:N:0:ATGAATCT',
		'TGGCGCACACAGTCGCTTTTAGACAA', '+', 'AAAAAEEEAEEEEEEE/EEEEEEEEE'
	],
	[
		'@NB501227:57:HGJM5BGX2:1:11101:12257:1061 1:N:0:ATGAAACT',
		'GTCATTTAGCTGTTCAATACATAGGG', '+', 'A/AAAE//EAA6EE//EAEEE/E//A'
	],
	[
		'@NB501227:57:HGJM5BGX2:1:11101:1717:1063 1:N:0:ATGAATCT',
		'CCACTGCCATTTCACTAATCGGAACC', '+', 'AAAAAAEEEEEEAEEEA/EEEEEEEE'
	],
	[
		'@NB501227:57:HGJM5BGX2:1:11101:13235:1063 1:N:0:ATTAATCT',
		'ACACCAATCAATCTCTCACACTATCC', '+', 'AAAAAE6EE/AEAEEEE/EEEEEEEE'
	]
];
is_deeply( [ map { $_->{'data'} } @{ @$result[0] } ],
	$exp, "right entries for reads R1" );

$exp = [
	[
		'@NB501227:57:HGJM5BGX2:1:11101:3408:1053 2:N:0:ATGAATCT',
'GGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN',
		'+',
'AA/###############################################################################################'
	],
	[
		'@NB501227:57:HGJM5BGX2:1:11101:4228:1060 2:N:0:ATGAATCT',
'ACCANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN',
		'+',
'AAA/##############################################################################################'
	],
	[
		'@NB501227:57:HGJM5BGX2:1:11101:12257:1061 2:N:0:ATGAAACT',
'GAAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN',
		'+',
'A///##############################################################################################'
	],
	[
		'@NB501227:57:HGJM5BGX2:1:11101:1717:1063 2:N:0:ATGAATCT',
'AACGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN',
		'+',
'AA//A#############################################################################################'
	],
	[
		'@NB501227:57:HGJM5BGX2:1:11101:13235:1063 2:N:0:ATTAATCT',
'AAACANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN',
		'+',
'6/A//#############################################################################################'
	]
];
is_deeply( [ map { $_->{'data'} } @{ @$result[1] } ],
	$exp, "right entries for reads R2" );

$exp = [
	[
		'@NB501227:57:HGJM5BGX2:1:11101:3408:1053 1:N:0:ATGAATCT', 'ATGAATCT',
		'+',                                                       'AAAAAEEE'
	],
	[
		'@NB501227:57:HGJM5BGX2:1:11101:4228:1060 1:N:0:ATGAATCT', 'ATGAATCT',
		'+',                                                       'AAAAA/EE'
	],
	[
		'@NB501227:57:HGJM5BGX2:1:11101:12257:1061 1:N:0:ATGAAACT',
		'ATGAAACT', '+', '//A///E/'
	],
	[
		'@NB501227:57:HGJM5BGX2:1:11101:1717:1063 1:N:0:ATGAATCT', 'ATGAATCT',
		'+',                                                       'AAAAAAEE'
	],
	[
		'@NB501227:57:HGJM5BGX2:1:11101:13235:1063 1:N:0:ATTAATCT',
		'ATTAATCT', '+', '////6666'
	]
];
is_deeply( [ map { $_->{'data'} } @{ @$result[2] } ],
	$exp, "right entries for reads I1" );

#print "\$exp = ".root->print_perl_var_def([map { $_->{'data'} } @{@$result[2]}] ).";\n";

#print "\$exp = ".root->print_perl_var_def([@values]).";\n";
