#!/bin/bash

#Required Software
module load samtools/1.20
module load bwa/0.7.17
module load pairtools/0.2.2

#Directory Structure
home=/home/FCAM/msmith
core=/core/projects/EBP/smith
scratch=/scratch/msmith
gid="intDF011"
IN_BAM="${scratch}/${gid}_alignedhic.bam"
CHROM_SIZES=${home}/yahs/contigs/intDF011.asm.hic.hap1.p_ctg.chrom.sizes
