#!/bin/bash

#Add this directory to your path so you can call the module scripts

##Assemble contigs
./hifiasm.sh
##Align Hi-C reads to contigs
./aln_hic.sh
##
