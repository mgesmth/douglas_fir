#!/bin/bash

if [[ ( $@ == "--help") ||  $@ == "-h" ]]
then
    echo ""
    echo "Usage: ./aln_hic.sh -a <THREADS_ALN> -s <THREADS_SORT -m <MEM> -r <REF.FA> -1 <READ_R1.fastq.gz> -2 <READ_R2.fastq.gz>"
    echo "                    -o <OUTDIR> -p <SAMPLENAME> -l <LIBRARYNAME>"
    echo ""
    echo "Align Hi-C reads to a genome assembly in preparation for contact map creation and/or manual curation."
    echo ""
    echo "Requirements:"
    echo "  bwa"
    echo "  samtools"
    echo ""
    echo "Please note script assumes read pairs following the naming convention _R1/_R2 and are gzipped fastqs."
    echo "If you have a different convention, either change the script or change your file names!"
    echo ""
    echo "-a <THREADS>             Threads for alignment (default 4)."
    echo "-s <THREADS>             Threads for sorting (default 4)."
    echo "-m <MEMORY>              Maximum memory per thread for sorting (default 768 MB)"
    echo "-r <REFERENCE.FA>        Path to assembly."
    echo "-1 <read_R1.fastq.gz>    Path to read 1."
    echo "-2 <read_R2.fastq.gz>    Path to read 2."
    echo "-o <OUTDIR>              Path to the output directory."
    echo "-p <SAMPLENAME>          A string containing the same name for RG tag. (Optional: default HiC_sample)"
    echo "-l <LIBRARYNAME>         A string containing the library name for RG tag. (Optional: default HiC_library)"
    echo ""
    echo ""
	exit 0
fi

#defaults
threads_aln=4
threads_sort=4
mem="768MB"
samp="HiC_sample"
lib="HiC_library"

OPTSTRING="a:s:m:1::2:o:p:l:r:"
while getopts ${OPTSTRING} opt
do
case ${opt} in
  a) threads_aln=${OPTARG};;
  s) threads_sort=${OPTARG};;
  m) mem=${OPTARG};;
  1) R1=${OPTARG};;
  2) R2=${OPTARG};;
  o) outdir=${OPTARG};;
  p) samp=${OPTARG};;
  l) lib=${OPTARG};;
  r) ref=${OPTARG}
  ?)
    echo "invalid option: -${opt}"
    exit 1 ;;
  esac
done

#required argument check
if [[ -z ${R1} || -z ${R2} || -z ${outdir} || -z ${ref} ]] ; then
  echo "[E]: Options -1, -2, -r and -o require arguments. Exiting."
  echo "[E]: Run ./aln_hic.sh -h or --help to see detailed usage."
  exit 1
fi

#Define additional variables
r1_string="_R1"
name=${R1//$r1_string/}
ref_name=$(basename ${ref})
outfix=$(echo "$name" | sed 's/fastq.gz/bam/')
rg="@RG\\tID:${name}\\tSM:${samp}\\tPL:LS454\\tLB:${lib}"

date
echo "[M]: Welcome. We are aligning ${R1} and ${R2} to ${ref_name}."

bwa mem -SP5M -t "$threads_aln" -R "$rg" "$ref" "$R1" "$R2" | \
samtools sort -n -@ "$threads_sort" -m "$mem" -O "bam" -o "${outdir}/${outfix}"

if [[ $? -eq 0 ]] ; then
  date
  echo "[M]: Alignment complete."
  exit 0
else
  echo "[E]: Alignment failed. Exit code $?"
  date
  exit 1
fi
