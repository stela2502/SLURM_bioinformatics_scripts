#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 1;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $bed_file, );
$bed_file ||= '';
my $exec = $plugin_path . "/../bin/polish_bed_like_files.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/polish_bed_like_files";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}


my $cmd =
    "perl -I $plugin_path/../lib  $exec.pl "
. " -bed_file " . $bed_file 
. " -drop_mt " # or not?
. " -debug";
#print "\$exp = ".root->print_perl_var_def($value ).";\n