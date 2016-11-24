#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 2;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $outfile, $infile, @options, );

my $exec = $plugin_path . "/../bin/reformat_bed_gtf_intersect.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/reformat_bed_gtf_intersect";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}

$infile = $plugin_path. "/data/test_new.xls";

@options = ( 'A', 'B');

my $cmd =
    "perl -I $plugin_path/../lib  $exec "
. " -outfile " . $outpath."/outfile"
. " -infile " . $infile 
. " -options " . join(' ', @options )
. " -names " # or not?
. " -debug";
system( $cmd );
#print "\$exp = ".root->print_perl_var_def($value ).";\n";