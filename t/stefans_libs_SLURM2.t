#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 3;
BEGIN { use_ok 'stefans_libs::SLURM' }

use FindBin;
my $plugin_path = $FindBin::Bin;

my ( $value, @values, $exp );
$value = { n => 2, N => 1, t => '0:00:01', A => 'snic2016-4-13' };
my $SLURM = stefans_libs::SLURM->new( $value );

is_deeply( ref($SLURM), 'stefans_libs::SLURM',
	'simple test of function stefans_libs::SLURM -> new()' );

my $module_str = $SLURM->load_SLURM_modules( 'GCC/4.9.3-2.25','OpenMPI/1.10.2', 'BEDTools/2.25.0', 'SAMtools/0.1.19');
print $module_str;
