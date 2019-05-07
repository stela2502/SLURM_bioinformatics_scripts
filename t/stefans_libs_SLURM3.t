#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 4;
BEGIN { use_ok 'stefans_libs::SLURM' }

use FindBin;
my $plugin_path = $FindBin::Bin;

my ( $value, @values, $exp );
#$value = { n => 2, N => 1, t => '0:00:01', A => 'snic2016-4-13' };
$value = {};
my $SLURM = stefans_libs::SLURM->new( $value );

is_deeply( ref($SLURM), 'stefans_libs::SLURM',
	'simple test of function stefans_libs::SLURM -> new()' );
## Here I want to check that one script that has no important lines in it should not be run.

is_deeply( -f "$plugin_path/data/Script_to_not_run.sh", 1, "file $plugin_path/data/Script_to_not_run.sh exists" );

is_deeply( $SLURM->runScript( "$plugin_path/data/Script_to_not_run.sh" ,'' ), -1, "Script is not run due to simplicity" );

#print "\$exp = ".root->print_perl_var_def($value ).";\n";
