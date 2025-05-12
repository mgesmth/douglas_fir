#!/bin/bash

#Assessing quality of both interior and coastal genome assemblies, and building contact maps

HIFI=/home/FCAM/msmith/hifi_data/intDF_allhifi.fastq.gz
HIC1=/home/FCAM/msmith/hiC_data/allhiC_R1.fastq.gz
HIC2=/home/FCAM/msmith/hiC_data/allhiC_R2.fastq.gz
INT_PRIM=/core/projects/EBP/smith/CBP_assemblyfiles/interior_primary_final.fa
INT_ALT=/core/projects/EBP/smith/CBP_assemblyfiles/interior_alternate_final.fa
COASTAL=/core/projects/EBP/smith/coastal/coastalDF_scaffrenamed_sorted.fa
scratch=/scratch/msmith

##Get kmer estimates of genome attributes (heterozygosity, repeat length, genome size, etc.) & kmer-based genome quality assessment
module load R/4.2.2 meryl/1.4.1 merqury/1.3
export PATH="/home/FCAM/msmith/R/x86_64-pc-linux-gnu-library/4.2:$PATH"
export PATH="/core/projects/EBP/smith/bin/genomescope2.0:$PATH"
kmers_out=/core/projects/EBP/smith/merqury_out/intDF_kmers
./modules/kmers.sh -t 12 -k 21 -h "$HIFI" -o "$kmers_out" -x -p "$INT_PRIM" -a "$INT_ALT"
module unload R/4.2.2 meryl/1.4.1 merqury/1.3

##Run Quast to get basic assembly stats (N50 etc.)
module load python/3.8.1 quast/5.2.0
quast_out=/home/FCAM/msmith/quast_out
./modules/quast.sh -t 12 -a "$INT_PRIM" -o "${quast_out}/int_prim"
./modules/quast.sh -t 12 -a "$INT_ALT" -o "${quast_out}/int_alt"
./modules/quast.sh -t 12 -a "$COASTAL" -o "${quast_out}/coastal"
module unload python/3.8.1 quast/5.2.0

#Hi-C contact maps - this requires a lot of resources, i.e. at least 36 threads and 1Tb RAM (as well as tmp disk space) to run reasonably quickly
#Note that RAM is set internally within these scripts and assumes at least 1Tb of RAM
module load bwa/0.7.17 samtools/1.20 pairtools/0.2.2 juicer/1.22.01 python/3.8.1
PREFIX="interior_primary_contacts"
./modules/alnhic_and_filter.sh "$INT_PRIM" "$HIC1" "$HIC2" "${scratch}/${PREFIX}" "$scratch"
PAIRS="${scratch}/${PREFIX}_nodups.pairs"
topdir=/core/projects/EBP/smith/juicer_primary
./modules/juicer_tools.sh -t 36 -d "$topdir" -c "$PAIRS" -g "$PREFIX" -z "$INT_PRIM" -o "${topdir}/${PREFIX}" -x "$scratch"
module unload bwa/0.7.17 samtools/1.20 pairtools/0.2.2 juicer/1.22.01 python/3.8.1
#From here, ran Juicebox on a local machine to control the setting of java variables. Code looked something like:
#java -Xms10G -Xms30G -jar Juicebox.jar "${PREFIX}.hic"

