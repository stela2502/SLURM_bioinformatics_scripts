#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 4;

use FindBin;
my $plugin_path = $FindBin::Bin;
my ( $value, $exp );
my $exec = $plugin_path . "/../bin/bowtie2_run_aurora.pl";
ok( -f $exec, "the script has been found" );

if ( -d "$plugin_path/data/bowtie2/" ) {
	system( 'rm ' . "$plugin_path/data/bowtie2/" );
}
if ( -f -f "$plugin_path/data/test_empty.fastq.sh" ) {
	unlink( -f "$plugin_path/data/test_empty.fastq.sh" );
}
$value =
    "perl -I $plugin_path/../lib $exec "
  . "-files $plugin_path/data/test_empty.fastq.gz "
  . "-options n 1 N 1 t '00:02:00' p 2 "
  . "-genome $plugin_path/data/hg38/hg38 "    ## does not even exist
  . "-coverage  $plugin_path/data/fake_hg38.chrom.sizes.txt "
  . "-bigwigTracks ~/nobackup/Jonas_Larsson/Roman_Galeev/BigWigTrackInfo.txt  -debug";
print "the command:\n$value\n";
system($value );

ok( -d "$plugin_path/data/bowtie2/", "the script created the outpath" );

ok( -f "$plugin_path/data/test_empty.fastq.sh",
	"the SBATCH script was created" );

open( IN, "<$plugin_path/data/test_empty.fastq.sh" );
$value = [ map { chomp; $_ } <IN> ];
close(IN);
$exp = [
	'#! /bin/bash',
	'#SBATCH -n 1',
	'#SBATCH -N 1',
	'#SBATCH -t 00:02:00',
	'#SBATCH -J test_empty.fastq',
	'#SBATCH -o test_empty.fastq%j.out',
	'#SBATCH -e test_empty.fastq%j.err',
"bowtie2 -x /home/slang/git_Projects/SLURM_bioinformatics_scripts/t/data/hg38/hg38  -p 2 -U '/home/slang/git_Projects/SLURM_bioinformatics_scripts/t/data/test_empty.fastq.gz' -S /home/slang/git_Projects/SLURM_bioinformatics_scripts/t/data/bowtie2/test_empty.fastq_bowtie2.sam",
'samtools view -Sb  /home/slang/git_Projects/SLURM_bioinformatics_scripts/t/databowtie2/test_empty.fastq.gz_bowtie2.sam | samtools sort -@ 1 - /home/slang/git_Projects/SLURM_bioinformatics_scripts/t/databowtie2/test_empty.fastq.gz_bowtie2.sorted',
'if  [ -f /home/slang/git_Projects/SLURM_bioinformatics_scripts/t/databowtie2/test_empty.fastq.gz_bowtie2.sorted.bam ]&&[ -s /home/slang/git_Projects/SLURM_bioinformatics_scripts/t/databowtie2/test_empty.fastq.gz_bowtie2.sorted.bam ]; then',
'rm -f /home/slang/git_Projects/SLURM_bioinformatics_scripts/t/databowtie2/test_empty.fastq.gz_bowtie2.sam',
	'fi',
	'## no -coverage option - no bigwig conversion',
	''
];

#print " \$exp = " . root->print_perl_var_def($value) . "\n ";

is_deeply( $value, $exp, "The script contains the right entries" )



