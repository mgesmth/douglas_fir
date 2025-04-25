#!/bin/bash

#Add this directory to your path so you can call the module scripts
#Most scripts have a help message, so call the script with -h or --help to see positional arguments

#RESOURCES:
####Max threads needed: 36
####Max RAM needed: 1Tb
####Time: >2 weeks (my server doesn't have time limits, so unsure about exact wall time).

HIFI=/home/FCAM/msmith/hifi_data/intDF_allhifi.fastq.gz
HIC_R1=/home/FCAM/msmith/hiC_data/allhiC_R1.fastq.gz
HIC_R2=/home/FCAM/msmith/hiC_data/allhiC_R2.fastq.gz
PREFIX="intDF011"
ARIMA_SITES="GATC,GANTC"

##Check quality of HiFi Reads

##Assemble contigs & convert primary GFA to FASTA
./hifiasm.sh "/home/FCAM/msmith/hifiasm_out/hifiasm1_1/${PREFIX}" "${HIC_R1}" "${HIC_R2}" "${HIFI}"
contig_prim_asm="/home/FCAM/msmith/hifiasm_out/hifiasm1_1/${PREFIX}.asm.hic.p_ctg.fasta"

##Align Hi-C reads to contigs and use pairtools to identify contacts, sort and dedup
./alnhic_andfilter.sh "${contig_prim_asm}" "${HIC_R1}" "${HIC_R2}" "${PREFIX}"
filtered_bam=/scratch/msmith/

##Scaffold contig assembly using validated Hi-C contacts
scaff_outdir=/core/projects/EBP/smith/scaffold
./scaffold.sh "${contig_prim_asm}" "${filtered_bam}" "${scaff_outdir}/${PREFIX}" "${ARIMA_SITES}"
scaff_prim_asm="${scaff_outdir}/${PREFIX}_scaffolds_final.fa"

##Run QUAST to get assembly quality statistics
quast_out=/home/FCAM/msmith/quast_out/1_3
./quast.sh "${scaff_prim_asm}" "${quast_out}"
