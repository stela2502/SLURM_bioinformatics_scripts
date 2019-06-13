#! /bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 00:02:00
#SBATCH -A lu2016-2-7
#SBATCH -p lu
#SBATCH -J testaf54
#SBATCH -o testaf54.%j.out
#SBATCH -e testaf54.%j.err
module purge
module load icc/2017.1.132-GCC-6.3.0-2.27 impi/2017.1.132 SAMtools/1.4.1 HISAT2/2.1.0 BEDTools/2.26.0 Java/1.8.0_92 
date
hisat2  -x /home/stefan/git/SLURM_bioinformatics_scripts/t/data/hg38/hg38 -U /home/stefan/git/SLURM_bioinformatics_scripts/t/data/test_empty.fastq.gz --threads 2 --add-chrname > $SNIC_TMP/test_empty.fastq_hisat.sam
date
samtools view -Sb  $SNIC_TMP/test_empty.fastq_hisat.sam | samtools sort -@ 1 -o $SNIC_TMP/test_empty.fastq_hisat.sorted.bam -
if  [ -f $SNIC_TMP/test_empty.fastq_hisat.sorted.bam ]&&[ -s $SNIC_TMP/test_empty.fastq_hisat.sorted.bam ]; then
rm -f $SNIC_TMP/test_empty.fastq_hisat.sam
fi
date
## here comes a subscript entry:
module purge
module load GCC/4.9.3-2.25 OpenMPI/1.10.2 icc/2016.1.150-GCC-4.9.3-2.25 impi/5.1.2.150 BEDTools/2.25.0 Java/1.8.0_72 picard/2.8.2 ucsc-tools/R2016a 
bedtools genomecov -bg -split -ibam $SNIC_TMP/test_empty.fastq_hisat.sorted.bam -g /home/stefan/git/SLURM_bioinformatics_scripts/t/data/fake_hg38.chrom.sizes.txt | sort -k1,1 -k2,2n > $SNIC_TMP/test_empty.fastq_hisat.bedGraph
polish_bed_like_files.pl -bed_file $SNIC_TMP/test_empty.fastq_hisat.bedGraph
bedGraphToBigWig $SNIC_TMP/test_empty.fastq_hisat.bedGraph /home/stefan/git/SLURM_bioinformatics_scripts/t/data/fake_hg38.chrom.sizes.txt $SNIC_TMP/test_empty.fastq_hisat.bw
if [ -f $SNIC_TMP/'test_empty.fastq_hisat.sorted.bam' ];then 
  cp $SNIC_TMP/test_empty.fastq_hisat.sorted.bam /home/stefan/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.sorted.bam
else
  ( >&2 echo 'file $SNIC_TMP/test_empty.fastq_hisat.sorted.bam was not produced - files in the same path:')
  for filename in $SNIC_TMP/*
  do
    (>&2 echo $filename) 
  done;
fi
if [ -f $SNIC_TMP/'test_empty.fastq_hisat.bw' ];then 
  cp $SNIC_TMP/test_empty.fastq_hisat.bw /home/stefan/git/SLURM_bioinformatics_scripts/t/data/HISAT2_OUT/test_empty.fastq_hisat.bw
else
  ( >&2 echo 'file $SNIC_TMP/test_empty.fastq_hisat.bw was not produced - files in the same path:')
  for filename in $SNIC_TMP/*
  do
    (>&2 echo $filename) 
  done;
fi
date
exit 0
