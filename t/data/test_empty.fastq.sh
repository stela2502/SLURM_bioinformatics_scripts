#! /bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 00:02:00
#SBATCH -J test_empty.fastq
#SBATCH -o test_empty.fastq%j.out
#SBATCH -e test_empty.fastq%j.err
#hisat2 -x /home/med-sal/git_Projects/SLURM_bioinformatics_scripts/t/data/hg38/hg38 -U /home/med-sal/git_Projects/SLURM_bioinformatics_scripts/t/data/test_empty.fastq.gz --threads 2 > /home/med-sal/git_Projects/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sam
#samtools view -Sb  /home/med-sal/git_Projects/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sam | samtools sort -@ 1 -o /home/med-sal/git_Projects/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sorted.bam -
#if  [ -f /home/med-sal/git_Projects/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sorted.bam ]&&[ -s /home/med-sal/git_Projects/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sorted.bam ]; then
#rm -f /home/med-sal/git_Projects/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sam
#fi
#bedtools genomecov -bg -split -ibam /home/med-sal/git_Projects/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sorted.bam -g /home/med-sal/git_Projects/SLURM_bioinformatics_scripts/t/data/fake_hg38.chrom.sizes.txt | sort -k1,1 -k2,2n > /home/med-sal/git_Projects/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.bedGraph
#bedGraphToBigWig is not working on aurora !?
#If it is working remove this warning!

