#! /bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 00:02:00
#SBATCH -A snic2016-4-13
#SBATCH -J test_empty.fastq
#SBATCH -o test_empty.fastq%j.out
#SBATCH -e test_empty.fastq%j.err
module load icc/2016.1.150-GCC-4.9.3-2.25 impi/5.1.2.150 HISAT2/2.0.4 BEDTools/2.25.0 Java/1.8.0_92 picard/2.8.2 BEDTools/2.25.0 
#hisat2 -x /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/hg38/hg38 -U /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/test_empty.fastq.gz --threads 2 --add-chrname > /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sam
#samtools view -Sb  /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sam | samtools sort -@ 1 -o /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sorted.bam -
#if  [ -f /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sorted.bam ]&&[ -s /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sorted.bam ]; then
#rm -f /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sam
#fi
#picard MarkDuplicates INPUT=/home/stefanl/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sorted.bam OUTPUT=/home/stefanl/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sorted_picard_deduplicated.bam REMOVE_DUPLICATES=TRUE M=/home/stefanl/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sorted_picard_metrix.txt
#bedtools genomecov -bg -split -ibam /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sorted_picard_deduplicated.bam -g /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/fake_hg38.chrom.sizes.txt | sort -k1,1 -k2,2n > /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sorted_picard_deduplicated.bedGraph

