#!/bin/bash

#Add this directory to your path so you can call the module scripts
#Most scripts have a help message, so call the script with -h or --help to see positional arguments

HIFI=/home/FCAM/msmith/hifi_data/intDF_allhifi.fastq.gz
HIC_R1=/home/FCAM/msmith/hiC_data/allhiC_R1.fastq.gz
HIC_R2=/home/FCAM/msmith/hiC_data/allhiC_R2.fastq.gz
PREFIX="intDF011"

##Assemble contigs & convert primary GFA to FASTA
./hifiasm.sh "/home/FCAM/msmith/hifiasm_out/hifiasm1_1/${PREFIX}" "${HIC_R1}" "${HIC_R2}" "${HIFI}"
contig_prim_asm="/home/FCAM/msmith/hifiasm_out/hifiasm1_1/${PREFIX}.asm.hic.p_ctg.fasta"
##Align Hi-C reads to contigs and use pairtools to identify contacts, sort and dedup
./alnhic_andfilter.sh "${contig_prim_asm}" "${HIC_R1}" "${HIC_R2}" "${PREFIX}"
chrom_sizes=
##Scaffold and Quast report
./pairtools.sh /home/FCAM/msmith/hifiasm_out/hifiasm1_1/intDF011.asm.bp.p_ctg.fasta intDF011
