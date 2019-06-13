#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 1;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $infile, $outfile, );

map { $_ ||= '' } $infile, $outfile;
my $exec = $plugin_path . "/../bin/dropUMI_from_bam.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/dropUMI_from_bam";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}


my $cmd =
    "perl -I $plugin_path/../lib  $exec.pl "
. " -infile " . $infile 
. " -outfile " . $outfile 
. " -debug";
#print "\$exp = ".root->print_perl_var_def($value ).";\n