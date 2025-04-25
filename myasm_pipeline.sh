#!/bin/bash

#Add this directory to your path so you can call the module scripts
#Most scripts have a help message, so call the script with -h or --help to see positional arguments

HIFI=/home/FCAM/msmith/hifi_data/intDF_allhifi.fastq.gz
HIC_R1=/home/FCAM/msmith/hiC_data/allhiC_R1.fastq.gz
HIC_R2=/home/FCAM/msmith/hiC_data/allhiC_R2.fastq.gz
PREFIX="intDF011"

##Assemble contigs
./hifiasm.sh "${prefix}" "${HIC_R1}" "${HIC_R2}" "${HIFI}"
##Align Hi-C reads to contigs
./aln_hic.sh /home/FCAM/msmith/hifiasm_out/hifiasm1_1/intDF011.asm.bp.p_ctg.fasta intDF011
##Use pairtools to identify valid ligation junctions, sort, and dedup Hi-C alignments
./pairtools.sh /home/FCAM/msmith/hifiasm_out/hifiasm1_1/intDF011.asm.bp.p_ctg.fasta intDF011
