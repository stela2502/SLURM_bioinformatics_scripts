#! /bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 00:02:00
#SBATCH -A snic2016-4-13
#SBATCH -p 2
#SBATCH -J test_empty.fastq
#SBATCH -o test_empty.fastq.%j.out
#SBATCH -e test_empty.fastq.%j.err
module purge
module load icc/2017.1.132-GCC-6.3.0-2.27 impi/2017.1.132 SAMtools/1.4.1 HISAT2/2.1.0 BEDTools/2.26.0 Java/1.8.0_92 
#hisat2  -x /projects/fs1/med-sal/git/SLURM_bioinformatics_scripts/t/data/hg38/hg38 -U /projects/fs1/med-sal/git/SLURM_bioinformatics_scripts/t/data/test_empty.fastq.gz --threads 2 --add-chrname > $SNIC_TMP/test_empty.fastq_hisat.sam
#samtools view -Sb  $SNIC_TMP/test_empty.fastq_hisat.sam | samtools sort -@ 1 -o $SNIC_TMP/test_empty.fastq_hisat.sorted.bam -
#if  [ -f $SNIC_TMP/test_empty.fastq_hisat.sorted.bam ]&&[ -s $SNIC_TMP/test_empty.fastq_hisat.sorted.bam ]; then
#rm -f $SNIC_TMP/test_empty.fastq_hisat.sam
#fi
#cp /projects/fs1/med-sal/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sorted.bam $SNIC_TMP/test_empty.fastq_hisat.sorted.bammodule purge
module load GCC/4.9.3-2.25 OpenMPI/1.10.2 icc/2016.1.150-GCC-4.9.3-2.25 impi/5.1.2.150 BEDTools/2.25.0 Java/1.8.0_72 picard/2.8.2 ucsc-tools/R2016a 
#picard MarkDuplicates INPUT=$SNIC_TMP/test_empty.fastq_hisat.sorted.bam OUTPUT=$SNIC_TMP/test_empty.fastq_hisat.sorted_picard_deduplicated.bam REMOVE_DUPLICATES=TRUE M=$SNIC_TMP/test_empty.fastq_hisat.sorted_picard_metrix.txt
#cp /projects/fs1/med-sal/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sorted_picard_deduplicated.bam $SNIC_TMP/test_empty.fastq_hisat.sorted_picard_deduplicated.bam
#
#bedtools genomecov -bg -split -ibam $SNIC_TMP/test_empty.fastq_hisat.sorted_picard_deduplicated.bam -g /projects/fs1/med-sal/git/SLURM_bioinformatics_scripts/t/data/fake_hg38.chrom.sizes.txt | sort -k1,1 -k2,2n > $SNIC_TMP/test_empty.fastq_hisat.sorted_picard_deduplicated.bedGraph
#cp /projects/fs1/med-sal/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sorted_picard_deduplicated.bedGraph $SNIC_TMP/test_empty.fastq_hisat.sorted_picard_deduplicated.bedGraph
#
#polish_bed_like_files.pl -bed_file $SNIC_TMP/test_empty.fastq_hisat.sorted_picard_deduplicated.bedGraph
#bedGraphToBigWig $SNIC_TMP/test_empty.fastq_hisat.sorted_picard_deduplicated.bedGraph /projects/fs1/med-sal/git/SLURM_bioinformatics_scripts/t/data/fake_hg38.chrom.sizes.txt $SNIC_TMP/test_empty.fastq_hisat.sorted_picard_deduplicated.bw
#cp /projects/fs1/med-sal/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sorted_picard_deduplicated.bw $SNIC_TMP/test_empty.fastq_hisat.sorted_picard_deduplicated.bw
#
#mv $SNIC_TMP/test_empty.fastq_hisat.sorted.bam /projects/fs1/med-sal/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sorted.bam
#mv $SNIC_TMP/test_empty.fastq_hisat.sorted_picard_deduplicated.bam /projects/fs1/med-sal/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sorted_picard_deduplicated.bam
#mv $SNIC_TMP/test_empty.fastq_hisat.sorted_picard_deduplicated.bw /projects/fs1/med-sal/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sorted_picard_deduplicated.bw
