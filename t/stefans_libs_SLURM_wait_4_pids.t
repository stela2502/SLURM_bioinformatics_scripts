#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 4;
BEGIN { use_ok 'stefans_libs::SLURM' }

use FindBin;
my $plugin_path = $FindBin::Bin;

my ( $value, @values, $exp );
$value = { n => 1, N => 1, t => '0:01:00', A => 'lsens2017-3-2' , 'p' => 'lsens' };

my $SLURM = stefans_libs::SLURM->new( $value );

is_deeply( ref($SLURM), 'stefans_libs::SLURM',
	'simple test of function stefans_libs::SLURM -> new()' );

my $cmd = "sleep 10"; ## does nothing just block a node for 10 sec.
my $outpath = "$plugin_path/data/output/SLURM_tmp";
unless ( -d $outpath){
	mkdir ( $outpath ) or die "I could not create the outpath\n$!\n";
}else {
	system( "rm -Rf $outpath/*" ); # get rid of ols scripts
}

my $start = time;

for ( my $i = 1; $i < 11; $i ++) {
	push(@values, $SLURM->run( $cmd, "$outpath/script_$i.sh" ) );
}

#print " \$exp = " . root->print_perl_var_def(\@values) . ";\n ";

is_deeply( scalar(@values), 10, "I got 10 ip's back?");

while( ! $SLURM->pids_finished (@values )) {
	print "sleep for 1 sec\n";
	sleep( 1 );
}

my $duration = time - $start;

ok ( $duration > 10, "more than 10 sec for the whole waiting: $duration sec\n" );

