#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 1;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $project_id, $access_token, $opath, );


map { $_ ||= '' } $project_id, $access_token, $opath;

my $exec = $plugin_path . "/../bin/download_fasta_from_Basespace_usingR.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/download_fasta_from_Basespace_usingR";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}


my $cmd =
    "perl -I $plugin_path/../lib  $exec.pl "
. " -project_id " . $project_id 
. " -access_token " . $access_token 
. " -opath " . $opath 
. " -debug";
#print "\$exp = ".root->print_perl_var_def($value ).";\n