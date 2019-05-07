#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 1;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $infile, $sequence, $outfile, );

map { $_ ||= ''} $infile, $sequence, $outfile;
my $exec = $plugin_path . "/../bin/filter_fastq_file.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/filter_fastq_file";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}


my $cmd =
    "perl -I $plugin_path/../lib  $exec.pl "
. " -infile " . $infile 
. " -sequence " . $sequence 
. " -outfile " . $outfile 
. " -debug";
#print "\$exp = ".root->print_perl_var_def($value ).";\n