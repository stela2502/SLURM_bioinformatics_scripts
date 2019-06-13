#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 2;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $data_path, @options, $path, $outpath, );

my $exec = $plugin_path . "/../bin/backupFastq.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/backupFastq";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}


my $cmd =
    "perl -I $plugin_path/../lib  $exec "
. " -data_path " . $data_path 
. " -options " . join(' ', @options )
. " -path " . $path 
. " -outpath " . $outpath 
. " -debug";
my $start = time;
system( $cmd );
my $duration = time - $start;
print "Execution time: $duration s\n";
#print "\$exp = ".root->print_perl_var_def($value ).";\n";