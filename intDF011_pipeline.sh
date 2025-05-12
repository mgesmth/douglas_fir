#!/bin/bash

#Add this directory to your path so you can call the module scripts
#Most scripts have a help message, so call the script with -h or --help to see positional arguments

#RESOURCES:
####threads needed: 36
####RAM needed: ~1Tb
####Time: >2 weeks (my server doesn't have time limits, so unsure about exact wall time).

scratch=/scratch/msmith
HIFI=/home/FCAM/msmith/hifi_data/intDF_allhifi.fastq.gz
HIC_R1=/home/FCAM/msmith/hiC_data/allhiC_R1.fastq.gz
HIC_R2=/home/FCAM/msmith/hiC_data/allhiC_R2.fastq.gz
PREFIX="intDF011"
ARIMA_SITES="GATC,GANTC"

##Check quality of HiFi Reads and Hi-C Reads
###Ran with threads=12 and mem=100G
module load fastqc/0.12.1 MultiQC/1.10.1 
fastqc_outdir_hifi=/home/FCAM/msmith/fastqc_out_hifi
fastqc_outdir_hic=/home/FCAM/msmith/fastqc_out_hic
./modules/fastqc.sh "${fastqc_outdir_hifi}" "${HIFI}"
./modules/fastqc.sh "${fastqc_outdir_hic}" "${HIC_R1}" "${HIC_R2}"
module unload fastqc/0.12.1 MultiQC/1.10.1 

##Assemble contigs & convert primary GFA to FASTA
###Ran with threads=36 and mem=800G
module load Hifiasm/0.20.0
./modules/hifiasm.sh "/home/FCAM/msmith/hifiasm_out/hifiasm1_1/${PREFIX}" "${HIC_R1}" "${HIC_R2}" "${HIFI}"
contig_prim_asm="/home/FCAM/msmith/hifiasm_out/hifiasm1_1/${PREFIX}.asm.hic.p_ctg.fasta"
module unload Hifiasm/0.20.0

##Align Hi-C reads to contigs and use pairtools to identify contacts, sort and dedup
###Ran with threads=36 and mem=1000G
module load samtools/1.20 bwa/0.7.17 pairtools/0.2.2
./modules/alnhic_andfilter.sh "${contig_prim_asm}" "${HIC_R1}" "${HIC_R2}" "${scratch}/${PREFIX}" "${scratch}"
filtered_bam=/scratch/msmith/"${PREFIX}_nodups.bam"
module unload samtools/1.20 bwa/0.7.17 pairtools/0.2.2

##Scaffold contig assembly using validated Hi-C contacts
###Ran with threads=36 and mem=700G
module load YaHS/1.2.2
scaff_outdir=/core/projects/EBP/smith/scaffold
./modules/scaffold.sh "${contig_prim_asm}" "${filtered_bam}" "${scaff_outdir}/${PREFIX}" "${ARIMA_SITES}"
scaff_prim_asm="${scaff_outdir}/${PREFIX}_scaffolds_final.fa"
module unload YaHS/1.2.2

##Run QUAST to get assembly quality statistics
###Ran with threads=12 and mem=150G
module load python/3.8.1 quast/5.2.0
quast_out=/home/FCAM/msmith/quast_out/1_3
./modules/quast.sh "${scaff_prim_asm}" "${quast_out}"
module unload python/3.8.1 quast/5.2.0

##Run Meryl, GenomeScope2 and Merqury to get kmer assembly quality statistics
module load R/4.2.2 meryl/1.4.1 merqury/1.3
export PATH="/home/FCAM/msmith/R/x86_64-pc-linux-gnu-library/4.2:$PATH"
kmer_out=/core/projects/EBP/smith/merqury_out
./modules/kmers.sh -t 18 -k 21 -h "${HIFI}" -o "${kmer_out}" -x -p "${scaff_prim_asm}"
module unload R/4.2.2 meryl/1.4.1 merqury/1.3
