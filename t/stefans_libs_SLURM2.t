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

my $module_str = $SLURM->load_SLURM_modules( 'GCC/4.9.3-2.25','OpenMPI/1.10.2', 'BEDTools/2.25.0', 'SAMtools/0.1.19');

$exp = "module load GCC/4.9.3-2.25 OpenMPI/1.10.2 BEDTools/2.25.0 SAMtools/0.1.19 ";
ok ($module_str eq  $exp, "module string created as expected '$module_str'\n");

@values = split( "\n",$SLURM->script( 'touch nothingatall.txt', 'useless' ));

$exp=["#! /bin/bash",
"#SBATCH -n 2",
"#SBATCH -N 1",
"#SBATCH -t 01:00:00",
"#SBATCH -A snic2016-4-13",
"#SBATCH -J useless",
"#SBATCH -o useless%j.out",
"#SBATCH -e useless%j.err",
"touch nothingatall.txt"];

is_deeply( \@values, $exp, "right script using the previousely saved options");

