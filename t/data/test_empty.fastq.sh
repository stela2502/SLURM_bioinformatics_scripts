#! /bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 00:02:00
#SBATCH -A snic2016-4-13
#SBATCH -J test_empty.fastq
#SBATCH -o test_empty.fastq%j.out
#SBATCH -e test_empty.fastq%j.err
#hisat2 -x /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/hg38/hg38 -U /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/test_empty.fastq.gz --threads 2 --add-chrname > /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sam
#samtools view -Sb  /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sam | samtools sort -@ 1 -o /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sorted.bam -
#if  [ -f /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sorted.bam ]&&[ -s /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sorted.bam ]; then
#rm -f /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sam
#fi
#bedtools genomecov -bg -split -ibam /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sorted.bam -g /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/fake_hg38.chrom.sizes.txt | sort -k1,1 -k2,2n > /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.bedGraph
#polish_bed_like_files.pl -bed_file /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.bedGraph
#bedGraphToBigWig /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.bedGraph /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/fake_hg38.chrom.sizes.txt /home/stefanl/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.bw

