#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 2;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, @options, $paired, $coverage, @files, $genome, $bigWigTracks, );

my $exec = $plugin_path . "/../bin/bowtie_run_aurora.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/bowtie_run_aurora";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}


my $cmd =
    "perl -I $plugin_path/../lib  $exec "
. " -options " . join(' ', @options )
. " -paired " . $paired 
. " -coverage " . $coverage 
. " -files " . join(' ', @files )
. " -genome " . $genome 
. " -bigWigTracks " . $bigWigTracks 
. " -debug";
my $start = time;
system( $cmd );
my $duration = time - $start;
print "Execution time: $duration s\n";
#print "\$exp = ".root->print_perl_var_def($value ).";\n";