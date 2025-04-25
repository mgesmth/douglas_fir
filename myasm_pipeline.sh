#!/bin/bash

#Add this directory to your path so you can call the module scripts

##Assemble contigs
./hifiasm.sh
##Align Hi-C reads to contigs
./aln_hic.sh /home/FCAM/msmith/hifiasm_out/hifiasm1_1/intDF011.asm.bp.p_ctg.fasta intDF011
##Use pairtools to identify valid ligation junctions, sort, and dedup Hi-C alignments
./pairtools.sh /home/FCAM/msmith/hifiasm_out/hifiasm1_1/intDF011.asm.bp.p_ctg.fasta intDF011
