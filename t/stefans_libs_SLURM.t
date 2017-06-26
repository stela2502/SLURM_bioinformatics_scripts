#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 7;
BEGIN { use_ok 'stefans_libs::SLURM' }

use FindBin;
my $plugin_path = $FindBin::Bin;

my ( $value, @values, $exp );
$value = { n => 2, N => 1, t => '01:00:00', A => 'snic2016-4-13' };
my $SLURM = stefans_libs::SLURM->new( $value );

is_deeply($SLURM->{'options'}->options(), { n => 2, N => 1, t => '01:00:00', A => 'snic2016-4-13' },"options stored correctly" );

ok( ref($SLURM) eq 'stefans_libs::SLURM',
	'simple test of function stefans_libs::SLURM -> new()' );

is_deeply( $value, {}, "SLURM variables removbed from hash");

$value = $SLURM->script( "some_command.pl", "my_name" );

#print "\$exp =  " . root->print_perl_var_def( [ split( "\n", $value ) ] ) . ";\n";
$exp = [
	'#! /bin/bash',
	'#SBATCH -n 2',
	'#SBATCH -N 1',
	'#SBATCH -t 01:00:00',
	'#SBATCH -A snic2016-4-13',
	'#SBATCH -J my_name',
	'#SBATCH -o my_name%j.out',
	'#SBATCH -e my_name%j.err',
	'some_command.pl'
];

is_deeply([ split( "\n", $value ) ], $exp, "right script created" );
$SLURM->{'debug'} = 1;
$value = root->filemap( $plugin_path."/data/output/stefans_libs_SLURM.a.b.c.sh" );
unlink( $value->{'total'} ) if ( -f  $value->{'total'} );
$SLURM-> run ( "some_command.pl", $value );

ok ( -f $value->{'total'}, 'script file created'. " ($value->{'total'})" );
open ( IN, "<$value->{'total'}" );
$value = [];
@values = <IN>;
close ( IN );
$value = [ map { chomp; $_ } @values ];
$exp = [
	'#! /bin/bash',
	'#SBATCH -n 2',
	'#SBATCH -N 1',
	'#SBATCH -t 01:00:00',
	'#SBATCH -A snic2016-4-13',
	'#SBATCH -J stefans_libs_SLURM.a.b.c',
	'#SBATCH -o stefans_libs_SLURM.a.b.c%j.out',
	'#SBATCH -e stefans_libs_SLURM.a.b.c%j.err',
	'some_command.pl'
];
is_deeply( $exp, $value, "the script file is as expected" );

$SLURM->{'options'}->save();

#print "\$exp = ".root->print_perl_var_def($value ).";\n";

