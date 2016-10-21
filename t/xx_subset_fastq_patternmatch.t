#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 2;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $infile, $outfile, $pattern, $where, );

my $exec = $plugin_path . "/../bin/subset_fastq_patternmatch.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/subset_fastq_patternmatch";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}


my $cmd =
    "perl -I $plugin_path/../lib  $exec.pl "
. " -infile " . $infile 
. " -outfile " . $outfile 
. " -pattern " . $pattern 
. " -where " . $where 
. " -debug";
#print "\$exp = ".root->print_perl_var_def($value ).";\n