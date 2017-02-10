#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 3;
use stefans_libs::flexible_data_structures::data_table;

use stefans_libs::SLURM;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, @options, @pids, $cmd, );

my $exec = $plugin_path . "/../bin/wait4pidsAndRun.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/wait4pidsAndRun";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}
mkdir ( $outpath );
@options = ( "sleep", '1' );

## first I need to create a SLURM process that I could monitor

my $SLURM= stefans_libs::SLURM->new( { 'A' => 'lu2016-2-7', 'n' => 1, 'N' => 1, 't'=> '00:00:20' } );

my $pid = $SLURM->run( 'sleep 10', "$outpath/sleep" );

ok( $pid != 1, "usable pid ($pid)");

sleep( 3 );
$cmd =
    "perl -I $plugin_path/../lib  $exec "
. " -options " . join(' ', @options )
. " -pids " . $pid
. " -cmd " . "'touch $outpath/awake.txt'" ;
system( $cmd );
sleep(1);
ok ( -f "$outpath/awake.txt", "downstream process output ($outpath/awake.txt)" );

#print "\$exp = ".root->print_perl_var_def($value ).";\n";