#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 2;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, @bed_files, $outfile, @options, );

my $exec = $plugin_path . "/../bin/SumUpBedFiles.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/SumUpBedFiles";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}


my $cmd =
    "perl -I $plugin_path/../lib  $exec "
. " -bed_files " . join(' ', @bed_files )
. " -outfile " . $outfile 
. " -options " . join(' ', @options )
. " -debug";
system( $cmd );
#print "\$exp = ".root->print_perl_var_def($value ).";\n";