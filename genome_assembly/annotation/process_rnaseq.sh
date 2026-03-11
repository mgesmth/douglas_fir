#!/bin/bash
#SBATCH -J fastqc
#SBATCH -p general
#SBATCH -q general
#SBATCH -c 24
#SBATCH --mem=150G
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err

set -e
date
echo "[M]: Host Name: `hostname`"

##Adapted from: https://gitlab.com/douglas-fir-transcriptome/De-novo-assembly-of-short-reads/-/blob/De-novo-assembly/

home=/home/FCAM/msmith
core=/core/projects/EBP/smith
scratch=/scratch/msmith
pe_dir=${core}/genome_annotation_shortread_data/pe_reads
raw_reads=${pe_dir}/raw
fastq_dir=${pe_dir}/fastqc_raw

module load fastqc/0.12.1 MultiQC/1.29

cd ${pe_dir}/raw
ls -1 *.fastq.gz > files.tmp
files=$(cat files.tmp)
echo -e "`date`:[M]: Generating FASTQC report on raw PE RNA-seq reads.\n" 
fastqc -o ${fastq_dir} -t 24 ${files}
cd ${fastq_dir}
multiqc .
rm ${pe_dir}/raw/files.tmp

echo -e "\n`date`:[M]: Finished generating FASTQC report of raw reads."
echo -e "`date`:[M]: Moving onto trimming adapters and removing low quality reads.\n"

outdir=${pe_dir}/trim
adaptors=${pe_dir}/NEBNext_dual_adaptors.fasta

module load Trimmomatic/0.39 java/22

cd ..
ls -1 *R1.fastq.gz > R1_files.tmp
for R1 in $(cat R1_files.tmp) ; do
  R2=$(echo "$R1" | sed 's/_R1/_R2/g')
  base=${R1/_R1.fastq.gz/}

  java -Xmx100G -jar $Trimmomatic PE \
  -threads 24 -phred33 -trimlog ${outdir}/${base}_log \
  ${R1} ${R2} \
  ${outdir}/${base}_trim_R1_paired.fastq.gz \
  ${outdir}/${base}_trim_R1_unpaired.fastq.gz \
  ${outdir}/${base}_trim_R2_paired.fastq.gz \
  ${outdir}/${base}_trim_R2_unpaired.fastq.gz \
  ILLUMINACLIP:${adaptors}:2:30:10:2:keepBothReads \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30
done && rm R1_files.tmp

echo -e "\n`date`:[M]: Done Trimmomatic."
echo -e "`date`:[M]: Generating FASTQC report for trimmed reads.\n"
fastq_dir=${pe_dir}/fastqc_trim

cd ${outdir}
ls -1 *_paired.fastq.gz > trim_files.tmp
files=$(cat trim_files.tmp)
fastqc -o ${fastq_dir} -t 24 ${files}
cd ${fastq_dir}
multiqc .
rm ${outdir}/trim_files.tmp

echo -e "\n`date`:[M]: Finished generating FASTQC report. Bye!"
