#! /usr/bin/perl -w

my $cmd = $ARGV[0];
my $name = $ARGV[1];
print join("\n",
'#! /bin/bash',
'#SBATCH -n 2',
'#SBATCH -N 1',
'#SBATCH -t 01:00:00',
'#SBATCH -J '.$name,
'#SBATCH -o '.$name.'_omp_%j.out',
'#SBATCH -e '.$name.'_omp_%j.err',
$cmd)

