#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 2;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, @options, $mapper_options, @files, $outpath, $fast_tmp, $genome, );

my $exec = $plugin_path . "/../bin/kallisto_run_aurora.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/kallisto_run_aurora";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}


my $cmd =
    "perl -I $plugin_path/../lib  $exec "
. " -options " . join(' ', @options )
. " -mapper_options " . $mapper_options 
. " -files " . join(' ', @files )
. " -outpath " . $outpath 
. " -fast_tmp " . $fast_tmp 
. " -genome " . $genome 
. " -debug";
my $start = time;
system( $cmd );
my $duration = time - $start;
print "Execution time: $duration s\n";
#print "\$exp = ".root->print_perl_var_def($value ).";\n";