#! /usr/bin/perl
use strict;
use warnings;
use Test::More tests => 26;
BEGIN { use_ok 'stefans_libs::file_readers::transfaq' }
use stefans_libs::file_readers::bed_file;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp );
my $OBJ = stefans_libs::file_readers::transfaq -> new({'debug' => 1});
is_deeply ( ref($OBJ) , 'stefans_libs::file_readers::transfaq', 'simple test of function stefans_libs::file_readers::transfaq -> new() ');

my $infile = "$plugin_path/data/test_zargos_out.txt";

ok ( -f $infile, "test file exists ($infile)" );

my $outpath = "$plugin_path/data/output/Zargos2pwm_plot";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}
mkdir($outpath); 

$OBJ->read_file ( $infile );

ok ($OBJ->{'filename'} eq $infile, "filename is strored" );

$OBJ->create_logo( $outpath, 0 );

ok ( -f "$outpath/1_1041.pdf", "logo R script" );
ok ( -f "$outpath/1_1041.pdf", "pdf logo" );

## to create them all in one go is tested in xx_Zargos2pwm_plot.t


my $peakfile = "$plugin_path/data/test_pirhanja_out.txt";
ok ( -f $peakfile, "test file 2 exists ($peakfile)" );
my $bed_file = stefans_libs::file_readers::bed_file->new();
$bed_file ->read_file ( $peakfile );

$value = $OBJ->map_motives_to_peaks ($bed_file,0 );
$value = $value -> sort ();
print $value->AsString();

#print "\$exp = ".root->print_perl_var_def($value ).";\n";


