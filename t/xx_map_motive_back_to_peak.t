#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 2;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $outfile, $peaks, $motives, );

my $exec = $plugin_path . "/../bin/iCLIP/map_motive_back_to_peak.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/map_motive_back_to_peak";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}
mkdir( $outpath );

$peaks = "$plugin_path/data/test_pirhanja_out.txt";
ok ( -f $peaks, 'peak file');
$motives = "$plugin_path/data/test_zargos_out.txt";
ok ( -f $motives, 'motives file');
$outfile = "$outpath/binding_sites.bedGraph";

my $cmd =
    "perl -I $plugin_path/../lib  $exec "
. " -outfile " . $outfile 
. " -peaks " . $peaks 
. " -motives " . $motives 
. " -debug";

system( $cmd );


ok( -f "$outpath/binding_sites.bedGraph.log", "log file created" );


#print "\$exp = ".root->print_perl_var_def($value ).";\n";