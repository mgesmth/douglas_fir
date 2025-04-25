#!/bin/bash

#Required Software
module load samtools/1.20
module load bwa/0.7.17

#Directory Structure
home=/home/FCAM/msmith
core=/core/projects/EBP/smith
scratch=/scratch/msmith
asm=$1
gid="intDF011"
hic1=${home}/hiC_data/allhiC_R1.fastq.gz
hic2=${home}/hiC_data/allhiC_R2.fastq.gz
#Output BAM is enormous, so please ensure you have enough disk (>1TB)
out="${scratch}/${gid}_alignedhic.bam"

samtools faidx ${asm}
cut -f1-2 "${asm}.fai" > "${asm}.chrom.sizes"
bwa index ${asm}
bwa mem -SP5M -t 36 "${asm}" "${hic1}" "${hic2}" | samtools view -bh -o "${out}"
