#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 2;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $adapter, @sample_barcode, $infile, $outfile, );

my $exec = $plugin_path . "/../bin/Filter_Bellodi_barcodes.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/Filter_Bellodi_barcodes";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}


my $cmd =
    "perl -I $plugin_path/../lib  $exec.pl "
. " -adapter " . $adapter 
. " -sample_barcode " . join(' ', @sample_barcode )
. " -infile " . $infile 
. " -outfile " . $outfile 
. " -debug";
#print "\$exp = ".root->print_perl_var_def($value ).";\n