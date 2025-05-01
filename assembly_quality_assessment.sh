#!/bin/bash

#Assessing quality of both interior and coastal genome assemblies, and building contact maps

HIFI=/home/FCAM/msmith/hifi_data/intDF_allhifi.fastq.gz
INT_PRIM=/core/projects/EBP/smith/CBP_assemblyfiles/interior_primary_final.fa
INT_ALT=/core/projects/EBP/smith/CBP_assemblyfiles/interior_alternate_final.fa
COASTAL=/core/projects/EBP/smith/coastal/coastalDF_scaffrenamed_sorted.fa

##Get kmer estimates of genome attributes (heterozygosity, repeat length, genome size, etc.) & kmer-based genome quality assessment
module load R/4.2.2 meryl/1.4.1 merqury/1.3
export PATH="/home/FCAM/msmith/R/x86_64-pc-linux-gnu-library/4.2:$PATH"
export PATH="/core/projects/EBP/smith/bin/genomescope2.0:$PATH"
kmers_out=/core/projects/EBP/smith/merqury_out/intDF_kmers
./kmers.sh -t 12 -k 21 -h "$HIFI" -o "$kmers_out" -x -p "$INT_PRIM" -a "$INT_ALT"
module unload R/4.2.2 meryl/1.4.1 merqury/1.3

##Run Quast to get basic assembly stats (N50 etc.)
module load python/3.8.1 quast/5.2.0
quast_out=/home/FCAM/msmith/quast_out
./quast.sh -t 12 -a "$INT_PRIM" -o "${quast_out}/int_prim"
./quast.sh -t 12 -a "$INT_ALT" -o "${quast_out}/int_alt"
./quast.sh -t 12 -a "$COASTAL" -o "${quast_out}/coastal"
module unload python/3.8.1 quast/5.2.0
