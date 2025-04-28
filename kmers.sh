#!/bin/bash

if [[ ( $@ == "--help") ||  $@ == "-h" ]]
then
    echo ""
    echo "Usage: ./kmers.sh -t <THREADS> -k <KMER_SIZE> -h <HIFI.fastq.gz> -o <OUTDIR> [-x -p <ASM1> -a <ASM2>]"
    echo ""
    echo "Get kmer statistics (including completeness and QV) for an assembly."
    echo ""
    echo "Requirements:"
    echo "  R > 4.0 with argparse, ggplot2, scales, minpack.lm"
    echo "	meryl"
    echo "  genomescope2.0"
    echo "  merqury"
    echo ""
    echo "PLEASE NOTE: I have hard-coded a path to a version of _submit_merqury.sh into this script, as I needed to" 
    echo "change it for my cluster. Please change/remove the path if using this script (variable is sub_merqury)."
    echo ""
    echo "-t <THREADS>         Threads (default 1)."
    echo "-k <KMER_SIZE>       Kmer size (default 21)."
    echo "-h <HIFI.fastq.gz>   Path to hifi reads."
    echo "-o <OUTPREFIX>       A string defining the output prefix (including path)."
    echo "-x                   Perform assembly assessment rather than just kmer profile of long reads (default: false)."
    echo "-p <ASM1.fa>         Path to primary/hap1 assembly (with -x)."
    echo "-a <ASM2.fa>         Path to alternate/hap2 assembly (with -x)." 
    echo ""
    echo ""
	exit 0
fi

asm_assessment="false"
threads=1

OPTSTRING="t:k:h:p:a:x:"
while getopts ${OPTSTRING} opt
do
case ${opt} in
  t) threads=${OPTARG};;
  k) k=${OPTARG};;
  h) hifi=${OPTARG};;
  o) output=${OPTARG};;
  p) prim=${OPTARG};;
  a) alt=${OPTARG};;
  x) asm_assessment="true" ;;
  ?)
    echo "invalid option: ${OPTARG}"
    exit 1 ;;
  esac
done

set errexit

#Check that required arguments are supplied
if [[ -z "${hifi}" || -z "${output}" ] ; then
  echo "[E]: Options -h and -o are required. Exiting."
  echo "[E]: Run ./kmers.sh -h or --help to see detailed usage."
  exit 1
fi

#Create meryl db for long reads
echo "[M]: Beginning kmer counting with meryl..."
meryl count k="${k}" threads="${threads}" "${hifi}" output "${output}.meryl"
echo "[M]: Done."
db="${output}".meryl

#Run genomescope2
echo "[M]: Running GenomeScope2.0..." 
meryl histogram threads="${threads}" k="${k}" "${db}" > "${output}.meryl.hist"
mkdir ${output}_gs2
gsoutdir="${output}_gs2"
genomescope.R -i "${output}.meryl.hist" -o "${gsoutdir}" -k "${k}"
echo "[M]: Done."

#Run assembly evaluation if selected (or exit if not)
if [[ -v "${prim}" && "${asm_assessment}" == "true" ]]; then
  sub_merqury="/home/FCAM/smith/_submit_merqury.sh"
  echo "[M]: Assembly evaluation mode has been selected."
  if [[ -v "${alt}" ]] ; then
    echo "[M]: Running assembly evaluation with two assemblies..."
    ${sub_merqury} "${output}.meryl" "${prim}" "${alt}" "${output}"
    echo "[M]: Done."
  else
    echo "[M]: Running assembly evaluation with one assembly..."
    ${sub_merqury} "${output}.meryl" "${prim}" "${output}"
    echo "[M]: Done."
  fi
  elif [[ -z "${prim}" && "${asm_assessment}" == "true" ]] ; then
    echo "[E]: Assembly evaluation mode has been triggered, but no primary assembly supplied. Exiting."
    echo "[E]: Run ./kmers.sh -h or --help to see detailed usage."
    exit 1
  elif [[ "${asm_assessment}" == "false" ]] ; then
    exit 0
fi
