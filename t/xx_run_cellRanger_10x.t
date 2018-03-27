#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 2;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, @modules, $indexFile, @options, $outpath, $tmp, );

my $exec = $plugin_path . "/../bin/run_cellRanger_10x.pl";
ok( -f $exec, 'the script has been found' );
$outpath = "$plugin_path/data/output/run_cellRanger_10x";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}

$indexFile = "$plugin_path/data/sample_seet.csv";
@modules = ("notOfInterest-NotrunAnyhow");
@options = qw(what count);

$tmp = "/tmp/cellRangerTMP/";


my $cmd =
    "perl -I $plugin_path/../lib  $exec "
. " -modules " . join(' ', @modules )
. " -indexFile " . $indexFile 
. " -options " . join(' ', @options )
. " -outpath " . $outpath 
. " -tmp " . $tmp 
. " -debug";
my $start = time;
system( $cmd );
my $duration = time - $start;
print "Execution time: $duration s\n";
#print "\$exp = ".root->print_perl_var_def($value ).";\n";