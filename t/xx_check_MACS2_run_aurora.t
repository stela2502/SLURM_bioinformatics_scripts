#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 4;

use FindBin;
my $plugin_path = $FindBin::Bin;
my ( $value, $exp, @tmp );
my $exec = $plugin_path . "/../bin/MACS2_run_aurora.pl";
ok( -f $exec, "the script has been found" );

if ( -d "$plugin_path/data/MACS2_out" ) {
	system("rm -Rf $plugin_path/data/MACS2_out ");
}

my $cmd =
    "perl -I $plugin_path/../lib  $exec "
  . "-files $plugin_path/data/test_empty.bam "
  . "-cfiles $plugin_path/data/test_empty.bam "
  . "-options f BAM g hs q 0.01 -debug";

system($cmd );

ok( -d "$plugin_path/data/MACS2_out", "outpath created" );
ok( -f "$plugin_path/data/MACS2_out/test_empty.rmdup.sh",
	"script file created" );

open( IN, "<$plugin_path/data/MACS2_out/test_empty.rmdup.sh" )
  or die "I could not open the new script file!\n";
@tmp = map { chomp; $_ } <IN>;

$exp = [
	'#! /bin/bash',
	'#SBATCH -n 1',
	'#SBATCH -N 1',
	'#SBATCH -t 02:00:00',
	'#SBATCH -J test_empty.rmdup',
	'#SBATCH -o test_empty.rmdup%j.out',
	'#SBATCH -e test_empty.rmdup%j.err',
"samtools rmdup $plugin_path/data/test_empty.bam $plugin_path/data/test_empty.rmdup.bam",
"samtools rmdup $plugin_path/data/test_empty.bam $plugin_path/data/test_empty.rmdup.bam",
"MACS2 -t $plugin_path/data/test_empty.rmdup.bam -c $plugin_path/data/test_empty.rmdup.bam "
	." -g hs -f BAM -q 0.01 -n $plugin_path/data/MACS2_out/test_empty.rmdup",
	''
  ]
;
#print " \$exp = " . root->print_perl_var_def( \@tmp ) . ";\n ";
is_deeply( \@tmp, $exp, "correct slurm script" );
