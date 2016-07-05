#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 5;

use FindBin;
my $plugin_path = $FindBin::Bin;
my ( $value, $exp, $cmd );
my $exec = $plugin_path . "/../bin/runCommand.pl";
ok( -f $exec, "the script has been found" );

if ( -d "$plugin_path/data/outpath/SLURM_run/" ) {
	system("rm -Rf $plugin_path/data/outpath/SLURM_run/");
}

$cmd =
" perl -I $plugin_path/../lib/ $exec -cmd 'Do nothing at all' -outfile '$plugin_path/data/outpath/SLURM_run/SLURM_result.txt' -I_have_loaded_all_modules -debug";

system($cmd );

ok( -f "$plugin_path/data/outpath/SLURM_run/SLURM_result.sh",
	"slurm script created" );

ok( -f "$plugin_path/data/outpath/SLURM_run/SLURM_result.txt.log",
	"perl log file created" );

open( F, "<$plugin_path/data/outpath/SLURM_run/SLURM_result.sh" )
  or die "I could not open the slurm script\n$!\n";
my @values = map { chomp; $_ } <F>;
close(F);

#print " \$exp = " . root->print_perl_var_def( \@values ) . ";\n ";

$exp = [
	'#! /bin/bash',
	'#SBATCH -n 10',
	'#SBATCH -N 1',
	'#SBATCH -t 02:00:00',
	'#SBATCH -J SLURM_result',
	'#SBATCH -o SLURM_result%j.out',
	'#SBATCH -e SLURM_result%j.err',
	'Do nothing at all'
];

is_deeply( \@values, $exp, "right entries in the SLURM script" );

$cmd .=
" -options n 10 N 2 'mail-user' someone\@somewhere.com 'mail-type' 'END'";

system($cmd );
open( F, "<$plugin_path/data/outpath/SLURM_run/SLURM_result.sh" )
  or die "I could not open the slurm script\n$!\n";
@values = map { chomp; $_ } <F>;
close(F);

#print " \$exp = " . root->print_perl_var_def( \@values ) . ";\n ";

#$exp = [
#	'#! /bin/bash',
#	'#SBATCH -n 10',
#	'#SBATCH -N 2',
#	'#SBATCH -t 02:00:00',
#	'#SBATCH -J SLURM_result',
#	'#SBATCH -o SLURM_result%j.out',
#	'#SBATCH -e SLURM_result%j.err',
#	'Do nothing at all'
#];

$exp = [
	'#! /bin/bash',
	'#SBATCH -n 10',
	'#SBATCH -N 2',
	'#SBATCH -t 02:00:00',
	'#SBATCH --mail-user someone@somewhere.com',
	'#SBATCH --mail-type END',
	'#SBATCH -J SLURM_result',
	'#SBATCH -o SLURM_result%j.out',
	'#SBATCH -e SLURM_result%j.err',
	'Do nothing at all'
];

is_deeply( \@values, $exp, "right entries in the SLURM script + mail option" );

#print " \$exp = " . root->print_perl_var_def( \@values ) . ";\n ";
