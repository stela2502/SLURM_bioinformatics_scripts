#! /bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 00:02:00
#SBATCH -J test_empty.fastq
#SBATCH -o test_empty.fastq%j.out
#SBATCH -e test_empty.fastq%j.err
bowtie2 -x /home/slang/git_Projects/SLURM_bioinformatics_scripts/t/data/hg38/hg38  -p 2 -U '/home/slang/git_Projects/SLURM_bioinformatics_scripts/t/data/test_empty.fastq.gz' -S /home/slang/git_Projects/SLURM_bioinformatics_scripts/t/data/bowtie2/test_empty.fastq_bowtie2.sam
samtools view -Sb  /home/slang/git_Projects/SLURM_bioinformatics_scripts/t/databowtie2/test_empty.fastq.gz_bowtie2.sam | samtools sort -@ 1 - /home/slang/git_Projects/SLURM_bioinformatics_scripts/t/databowtie2/test_empty.fastq.gz_bowtie2.sorted
if  [ -f /home/slang/git_Projects/SLURM_bioinformatics_scripts/t/databowtie2/test_empty.fastq.gz_bowtie2.sorted.bam ]&&[ -s /home/slang/git_Projects/SLURM_bioinformatics_scripts/t/databowtie2/test_empty.fastq.gz_bowtie2.sorted.bam ]; then
rm -f /home/slang/git_Projects/SLURM_bioinformatics_scripts/t/databowtie2/test_empty.fastq.gz_bowtie2.sam
fi
## no -coverage option - no bigwig conversion

