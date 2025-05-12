#!/bin/bash


if [[ ( $@ == "--help") ||  $@ == "-h" ]]
then
    echo "Usage: ./juicer_tools.sh -t <THREADS> -d <TOPDIR> -c <PAIRS> -g <GENOMEID> -z <GENOMEPATH> -o <OUTPUT> [-p <CHROMSIZES> -x <TMPDIR>]"
    echo ""
    echo "Build contact maps from aligned Hi-C data."
    echo ""
    echo "dependencies:"
    echo ""
    echo "    java"
    echo "    samtools"
    echo "    python3"
    echo "    juicer"
    echo ""
    echo "Requires at least 4 cores and 64GB RAM."
    echo "Note this script is not currently built for fragment-level hic contact matrices."
    echo ""
    echo "positional arguments:"
    echo ""
    echo "-t <THREADS>     Number of threads."
    echo "-d <TOPDIR>      Top level directory. Must contain ./fastq, which contains Hi-C fastqs or is soft-linked to them and ./references, with the reference (or soft-linked to it)."
    echo "-c <PAIRS>       Hi-C contacts in pairs format (outputted by pairtools)."
    echo "-g <GENOMEID>    A unique identifier for your genome."
    echo "-z <GENOMEPATH>  Path to the reference genome."
    echo "-p <CHROMSIZES>  Path to chromosome size file (optional; can be produced internally)."    
    echo "-o <OUTPUT>      Prefix for output files."
    echo "-x <TMPDIR>      Path to temporary directory for temp files (optional)."
    echo ""
    echo ""
	exit 0
fi

#Defaults
chromsizes=""
tmpdir="."
threads=4

OPTSTRING="t:d:c:g:p:z:o:x:"
while getopts ${OPTSTRING} opt
do
    case ${opt} in
	t)
	 threads=${OPTARG};;
	d)
	 topdir=${OPTARG};;
	c)
	 contacts=${OPTARG};;
	g)
	 genome_id=${OPTARG};;
	p)
	 chromsizes=${OPTARG};;
	z)
	 ref=${OPTARG};;
	o)
	 output=${OPTARG};;
	x)
	 tmpdir=${OPTARG};;
  :)
    echo "option -${opt} requires an argument."
    exit 1
	;;
  ?)
    echo "invalid option: ${opt}"
    exit 1
	;;
    esac
done

#Check required options
if [[ -z "$topdir" || -z "$contacts" || -z "$genome_id" || -z "$ref" || -z "$output" ]] ; then
  echo "[E]: Options -d,-c,-g,-z, and -o require arguments. Exiting 1."
  echo "[E]: Run ./juicer_tools.sh -h or --help for detailed usage."
  exit 1
fi

set -e

#Chromsizes file
cd ${topdir}/references

if [ ! -f "${assembly}.fai" ]; then
	echo "[M]: Faidx not found. Indexing..."
	samtools faidx ${assembly}
	echo "[M]: Done."
	echo "[M]: Creating chrom.sizes file..."
	cut -f1-2 "${assembly}.fai" > "${asm_prefix}.chrom.sizes"
	echo "[M]: Done."
  rm "${assembly}.fai"
else
	echo "[M]: Faidx found."
	if [ ! -f "${asm_prefix}.chrom.sizes" ]; then
		echo "[M]: Chrom.sizes file not found. Creating..."
		cut -f1-2 "${assembly}.fai" > "${asm_prefix}.chrom.sizes"
		echo "[M]: Done."
	else
		echo "[M]: Chrom.sizes file found."
	fi
fi

CHROM_SIZES="${asm_prefix}.chrom.sizes"

#this horrible section is to check if the .pairs file is correctly formatted with dummy fragment values and tab field seps.
fieldcheck=`awk '/^#/ {print $0} !/^#/ {exit}' "${contacts}" | grep "columns" | grep "frag1"`
if [ -z "$fieldcheck" ] ; then
	if [ -f "${tmpdir}/contacts_corrected.pairs" ]; then
	    echo "[M]: Original .pairs file not formatted correctly, but corrected temp file was found. Continuing with this file."
	    contacts="${tmpdir}/contacts_corrected.pairs"
	else
      echo "[M]: Beginning reformatting of .pairs file..."
	    awk '
	    BEGIN { OFS = "\t" }
	    /^#/ {
    	    if ($0 ~ /^#columns:/) {
        	sub(/pair_type/, "frag1\tfrag2");
        	print;
    	    } else {
        	print;
    	    }
    		next;
	    }	
	    {
    		# Convert + to 0 and - to 1 in fields 6 and 7
    		if ($6 == "+") $6 = 0;
    		else if ($6 == "-") $6 = 1;

    		if ($7 == "+") $7 = 0;
    		else if ($7 == "-") $7 = 1;

    		# Build a list of fields from $1 to $(NF-1)
    		out = $1;
    		for (i = 2; i < NF; i++) {
        	out = out OFS $i;
    	    }
    	    print out, 0, 1;  # This uses OFS="\t" correctly
	    }' "${contacts}" > "${tmpdir}/contacts_corrected.pairs"
	    #reset contacts
	    contacts="${tmpdir}/contacts_corrected.pairs"
	    echo "[M]: Done."
	fi
else
	echo "[M]: Pairs file assumed correctly formatted."
fi

#Finally, make the .hic file
echo "[M]: Beginning .hic file creation."
java -XX:+UseParallelGC -Xms250G -Xmx400G -jar $JUICER pre -v --threads "${threads}" -t "${tmpdir}" "${contacts}" "${output}.hic" "${chromsizes}"
if [ $? -eq 0 ]; then
  echo "[M]: juicer_tools pre succeeded."
  exit 0
else
  echo "[E]: juicer_tools pre failed. Exiting."
  exit 1
fi 
