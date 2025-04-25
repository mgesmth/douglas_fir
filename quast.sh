#!/bin/bash

if [[ ( $@ == "--help") ||  $@ == "-h" ]]
then
    echo ""
    echo "Usage: ./quast.sh <ASSEMBLY> <OUT_PREFIX>"
    echo ""
    echo "Requirements:"
    echo "	bwa"
    echo "	samtools"
    echo "	pairtools"
    echo ""
    echo "<ASSEMBLY>        Path to assembly to be assessed."
    echo "<OUT_PREFIX>      A prefix for output, including path to out directory."
    echo ""
    echo "I recommend sending everything to scratch as these files can be massive."
    echo ""
	exit 0
fi

#Required Software
module load quast/5.2.0

asm=$1
out=$2

quast.py -t 12 --split-scaffolds --large -o ${out} ${asm}
