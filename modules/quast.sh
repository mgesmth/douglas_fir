#!/bin/bash

if [[ ( $@ == "--help") ||  $@ == "-h" ]]
then
    echo ""
    echo "Usage: ./quast.sh -t <THREADS> -a <ASSEMBLY> -o <OUTDIR>"
    echo ""
    echo "Requirements:"
    echo "	Python > 3"
    echo "	Quast"
    echo ""
    echo "-t <THREADS>	    Threads. Default 1."
    echo "-a <ASSEMBLY>     Path to assembly to be assessed."
    echo "-o <OUTDIR>       Output directory."
    echo ""
    echo ""
	exit 0
fi

OPTSTRING="t:a:o:"
while getopts ${OPTSTRING} opt
do
case ${opt} in
  t) threads=${OPTARG};;
  a) asm=${OPTARG};;
  o) outdir=${OPTARG};;
  ?)
    echo "invalid option: ${OPTARG}"
    exit 1 ;;
  esac
done

python3 quast.py -t ${threads} --split-scaffolds --large -o ${outdir} ${asm}
