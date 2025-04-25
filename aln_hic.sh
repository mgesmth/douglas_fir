#!/bin/bash

if [[ ( $@ == "--help") ||  $@ == "-h" ]]
then
    echo "USAGE: ./aln_hic.sh <ASM.fa> <HIC_R1.fastq.gz> <HIC_R2.fastq.gz> <GENOMEID>"
    echo ""
    echo "Align Hi-C reads to a genome assembly."
    echo ""
    echo "Positional arguments:"
    echo "<ASM>                The path to the assembly (FASTA)."
    echo "<HIC_R1.fastq.gz>    The path to Hi-C reads R1."
    echo "<HIC_R2.fastq.gz>    The path to Hi-C reads R2."
    echo "<GENOMEID>           A unique string to identify this assembly."
    echo ""
    echo ""
     exit 0
fi
    
#Required Software
module load samtools/1.20
module load bwa/0.7.17

#Directory Structure
home=/home/FCAM/msmith
core=/core/projects/EBP/smith
scratch=/scratch/msmith
asm=$1
gid=$2
hic1=$3
hic2=$4
#Output BAM is enormous, so please ensure you have enough disk (>1TB)
out="${scratch}/${gid}_alignedhic.bam"

samtools faidx ${asm}
cut -f1-2 "${asm}.fai" > "${asm}.chrom.sizes"
bwa index ${asm}
bwa mem -SP5M -t 36 "${asm}" "${hic1}" "${hic2}" | samtools view -bh -o "${out}"
