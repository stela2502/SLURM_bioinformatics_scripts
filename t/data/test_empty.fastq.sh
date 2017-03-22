#! /bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 00:02:00
#SBATCH -A snic2016-4-13
#SBATCH -J test_empty.fastq
#SBATCH -o test_empty.fastq%j.out
#SBATCH -e test_empty.fastq%j.err
module purge
module load icc/2015.3.187-GNU-4.9.3-2.25 impi/5.0.3.048 Bowtie/1.1.2 
#bowtie -m 500 -v2 -best -strata -p 2 --sam /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/hg38/hg38 '/home/stefanl/git/SLURM_bioinformatics_scripts/t/data/test_empty.fastq.gz' > /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/bowtie/test_empty.fastq_bowtie.sam
module purge
module load icc/2016.1.150-GCC-4.9.3-2.25 impi/5.1.2.150 SAMtools/1.3.1 BEDTools/2.25.0 ucsc-tools/R2016a 
#samtools view -Sb  /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/bowtie/test_empty.fastq_bowtie.sam | samtools sort -@ 1 -o /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/bowtie/test_empty.fastq_bowtie.sorted.bam -
#if  [ -f /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/bowtie/test_empty.fastq_bowtie.sorted.bam ]&&[ -s /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/bowtie/test_empty.fastq_bowtie.sorted.bam ]; then
#rm -f /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/bowtie/test_empty.fastq_bowtie.sam
#fi

#bedtools genomecov -bg -split -ibam /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/bowtie/test_empty.fastq_bowtie.sorted.bam -g /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/fake_hg38.chrom.sizes.txt | sort -k1,1 -k2,2n > /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/bowtie/test_empty.fastq_bowtie.bedGraph

#polish_bed_like_files.pl -bed_file /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/bowtie/test_empty.fastq_bowtie.bedGraph
#bedGraphToBigWig /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/bowtie/test_empty.fastq_bowtie.bedGraph /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/fake_hg38.chrom.sizes.txt /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/bowtie/test_empty.fastq_bowtie.bw

