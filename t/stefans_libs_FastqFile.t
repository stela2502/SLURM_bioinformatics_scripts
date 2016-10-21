#! /usr/bin/perl
use strict;
use warnings;
use Test::More tests => 9;
BEGIN { use_ok 'stefans_libs::FastqFile' }
use stefans_libs::root;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp );
my $OBJ = stefans_libs::FastqFile -> new({'debug' => 1});
is_deeply ( ref($OBJ) , 'stefans_libs::FastqFile', 'simple test of function stefans_libs::FastqFile -> new() ');

#print "\$exp = ".root->print_perl_var_def($value ).";\n";
$value = $OBJ -> open_file ($plugin_path."/data/repeat.fastq.gz" );
is_deeply ( ref($value), 'GLOB', 'opened gzipped txt file' );
@values = <$value>;
close ( $value );
#print "\$exp = ".root->print_perl_var_def([@values[0..3]]).";\n";
$exp = [ '@M04223:22:000000000-AUY26:1:1101:15371:1335 1:N:0:1
', 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
', '+
', '1>>>>>>1>>>>AAAAAA>>>>>>>><<<<B@@<<<<<<<:::99999999
' ];
is_deeply ( [@values[0..3]], $exp, "the right sequnce has been read\n" );


$value = $OBJ -> open_file ($plugin_path."/data/repeat.fastq" );

is_deeply ( ref($value), 'GLOB', 'opened normal txt file' );
@values = <$value>;
close ( $value );
is_deeply ( [@values[0..3]], $exp, "the right sequnce has been read\n" );


$OBJ -> exclude_read ($plugin_path."/data/repeat.fastq.gz", 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT', $plugin_path."/data/outpath/repeat.fastq" );
$value = $OBJ -> open_file ($plugin_path."/data/outpath/repeat.fastq" );
@values = <$value>;
close ( $value );
#print "\$exp = ".root->print_perl_var_def([@values]).";\n";
$exp = [ '@M04223:22:000000000-AUY26:1:1101:16818:1622 1:N:0:1
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
' ];

is_deeply ( [@values], $exp, "the right sequnce has been read\n" );


## and now test the extract_barcode function

#TGTGGCACGATGTGACGAGGCGGCCGAGTGGTTAAGGCGATGGACTGCTAA

my $adapter =  'AGATCGGAAGAGCGGTTCAG';
my $sample_barcode = 'GGCA';


$OBJ -> extract_barcode ( $plugin_path."/data/repeat.fastq.gz",  $adapter, $sample_barcode, 0, $plugin_path."/data/outpath/Barcodes_GGCA.fastq"  );
$value = $OBJ -> open_file ( $plugin_path."/data/outpath/Barcodes_GGCA.fastq" );
@values = <$value>;
close ( $value );
#print "\$exp = ".root->print_perl_var_def([@values]).";\n";
$exp = [ '@M04223:22:000000000-AUY26:1:1101:15584:1697:TGTCG 1:N:0:1
', 'ATGTGACGAGGCGGCCGAGTGGTTAAGGCGATGGCGGTTCAG
', '+
', '1>1F13A1100A00AA///F//FFB111A//E//0BGAEGBB
', '@M04223:22:000000000-AUY26:1:1101:12954:1698:AD11:TTCAG 1:N:0:1
', 'TCAAGCAGTTAACGCGGGTTCGATTCCCGGG
', '+
', 'BB@311A1FGEDE0A0?CEGF?AGHFDF???
' ];

is_deeply(\@values, $exp, "right sequences returned" );

$OBJ -> select_4_str ( $plugin_path."/data/outpath/Barcodes_GGCA.fastq", ":AD\\d",0, $plugin_path."/data/outpath/Barcodes_GGCA_AD.fastq" );
$value = $OBJ -> open_file ( $plugin_path."/data/outpath/Barcodes_GGCA_AD.fastq" );
@values = <$value>;
close ( $value );
$exp = [ '@M04223:22:000000000-AUY26:1:1101:12954:1698:AD11:TTCAG 1:N:0:1
', 'TCAAGCAGTTAACGCGGGTTCGATTCCCGGG
', '+
', 'BB@311A1FGEDE0A0?CEGF?AGHFDF???
' ];
is_deeply(\@values, $exp, "right select_4_str sequences returned" );



