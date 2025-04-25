#!/bin/bash

if [[ ( $@ == "--help") ||  $@ == "-h" ]]
then
    echo ""
    echo "Usage: ./scaffold.sh <ASSEMBLY> <FILT_CONTACTS.BAM> <OUT_PREFIX> <SITES>"
    echo ""
    echo "Requirements:"
    echo "	YaHS"
    echo ""
    echo "<ASSEMBLY>            Path to contig assembly."
    echo "<FILT_CONTACTS.BAM>   Path to filtered contacts in BAM (.pa5 also works)."
    echo "<OUT_PREFIX>          A prefix for output, including path to out directory."
    echo "<SITES>               A string of the restriction sites used in Hi-C library prep."
    echo ""
    echo ""
	exit 0
fi

#Required Software
module load YaHS/1.2.2

asm=$1
bam=$2
out=$3
sites=$4

yahs --no-contig-ec -l 10 -e "${sites}" -o ${out} ${asm} ${bam}

