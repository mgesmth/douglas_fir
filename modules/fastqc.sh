#!/bin/bash

if [[ ( $@ == "--help") ||  $@ == "-h" ]]
then
    echo ""
    echo "Usage: ./fastqc.sh <OUTDIR> <READS1> <READS2>"
    echo ""
    echo "Requirements:"
    echo "	FastQC"
    echo ""
    echo "<OUTDIR>        Path to out directory."
    echo "<READS1>        First reads file to be evaluated."
    echo "<READS2>        Second reads file to be evaluated (optional)."
    echo ""
    echo ""
	exit 0
fi

if [ $# -eq 2 ] ; then
  outdir=$1
  read=$2
  cd ${outdir}
  fastqc -o ./ -t 12 -d ${read1}
  multiqc .
elif [ $# -eq 3 ] ; then
  outdir=$1
  read1=$2
  read2=$3
  cd ${outdir}
  fastqc -o ./ -t 12 -d ${read1} ${read2}
  multiqc .
else
  echo "Unknown number of parameters. Exiting."
  exit 1
fi
  
