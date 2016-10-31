#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 2;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $outpath, @options, $infile, );

$infile = "$plugin_path/data/test_zargos_out.txt";

ok ( -f $infile, "test file exists ($infile)" );

$outpath = "$plugin_path/data/output/Zargos2pic";

my $exec = $plugin_path . "/../bin/Zargos2pwm_plot.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/Zargos2pwm_plot";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}

@options= qw(sampleName data1);

my $cmd =
    "perl -I $plugin_path/../lib  $exec "
. " -outpath " . $outpath 
. " -options " . join(' ', @options )
. " -infile " . $infile 
. " -debug";

system( $cmd );

ok ( -d $outpath, "outpath created ($outpath)" );

#print "\$exp = ".root->print_perl_var_def($value ).";\n";