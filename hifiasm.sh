#!/bin/bash

hifidir=/home/FCAM/msmith/hifi_data
hiCdir=/home/FCAM/msmith/hiC_data
outdir=/home/FCAM/msmith/hifiasm_out/hifiasm1_1

module load Hifiasm/0.20.0

hifiasm -o $outdir/intDF011.asm -t 36 \
--h1 $hiCdir/allhiC_R1.fastq.gz --h2 $hiCdir/allhiC_R2.fastq.gz \
-f 39 \
$hifidir/intDF_allhifi.fastq.gz

#Transform resulting primary GFA to fasta
awk '/^S/{print ">"$2"\n"$3}' ${outdir}/intDF011.asm.bp.p_ctg.gfa | fold > ${outdir}/intDF011.asm.bp.p_ctg.fasta
#in GFA files, the "segment" (i.e. actual sequence, or contig) starts with S. So finding all lines starting with
#S, print >, col $2 (which is the contig name) and then on the next line, $3 (the sequence)
#fold line wraps the file
