#!/bin/bash

hifidir=/home/FCAM/msmith/hifi_data
hiCdir=/home/FCAM/msmith/hiC_data
outdir=/home/FCAM/msmith/hifiasm_out/hifiasm1_1

module load Hifiasm/0.20.0

hifiasm -o $outdir/intDF010.asm -t 36 \
--h1 $hiCdir/allhiC_R1.fastq.gz --h2 $hiCdir/allhiC_R2.fastq.gz \
-f 39 \
$hifidir/intDF_allhifi.fastq.gz
