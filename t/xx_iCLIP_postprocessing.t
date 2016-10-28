#! /usr/bin/perl
system( "module load GCCcore/4.9.3 GCC/4.9.3-2.25 GSL/2.1 ucsc-tools/R2016a" );
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 1;
use stefans_libs::flexible_data_structures::data_table;


use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, @options, @infiles, $outpath, );

my $exec = $plugin_path . "/../bin/iCLIP_postprocessing.pl";
ok( -f $exec, 'the script has been found' );
$outpath = "$plugin_path/data/output/iCLIP_postprocessing";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}


my $cmd =
    "perl -I $plugin_path/../lib  $exec.pl "
. " -options " . join(' ', @options )
. " -infiles " . join(' ', @infiles )
. " -outpath " . $outpath 
. " -debug";
#print "\$exp = ".root->print_perl_var_def($value ).";\n