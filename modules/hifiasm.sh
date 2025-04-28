#!/bin/bash

if [[ ( $@ == "--help") ||  $@ == "-h" ]]
then
    echo "USAGE: ./hifiasm.sh <GENOMEID> <HIC_R1.fastq.gz> <HIC_R2.fastq.gz> <HIFIREADS.fastq.gz>"
    echo ""
    echo "Use hifiasm integrated mode to assemble contigs."
    echo ""
    echo "Positional arguments:"
    echo "<GENOMEID>             A unique string to identify this assembly, prefixed with the path to desired output directory."
    echo "<HIC_R1.fastq.gz>      The path to Hi-C reads R1."
    echo "<HIC_R2.fastq.gz>      The path to Hi-C reads R2."
    echo "<HIFIREADS.fastq.gz>   The path to PacBio HiFi long reads (all files merged into one)."
    echo ""
    echo ""
     exit 0
fi

#Required Software
module load Hifiasm/0.20.0

gid="$1"
hic1=$2
hic2=$3
hifi=$4


hifiasm -o ${gid}.asm -t 36 \
--h1 ${hic1} --h2 ${hic2} \
-f 39 \
${hifi}

#Transform resulting primary GFA to fasta
awk '/^S/{print ">"$2"\n"$3}' ${outdir}/intDF011.asm.bp.p_ctg.gfa | fold > ${outdir}/intDF011.asm.bp.p_ctg.fasta
#in GFA files, the "segment" (i.e. actual sequence, or contig) starts with S. So finding all lines starting with
#S, print >, col $2 (which is the contig name) and then on the next line, $3 (the sequence)
#fold line wraps the file
