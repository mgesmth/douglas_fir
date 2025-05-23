#!/bin/bash

if [[ ( $@ == "--help") ||  $@ == "-h" ]]
then
    echo ""
    echo "Usage: ./alnhic_andfilter.sh <ASSEMBLY> <HIC1> <HIC2> <OUT_PREFIX> <TMPDIR>"
    echo ""
    echo "Requirements:"
    echo "	bwa"
    echo "	samtools"
    echo "	pairtools"
    echo ""
    echo "<ASSEMBLY>        Path to assembly."
    echo "<HIC1>            Path to Hi-C R1 in .fastq.gz format."
    echo "<HIC2>            Path to Hi-C R2 in .fastq.gz format."
    echo "<OUT_PREFIX>      A prefix for output, including path to out directory."
    echo "<TMPDIR>          A temporary directory to send sorting files to, i.e. scratch."
    echo ""
    echo "I recommend sending everything to scratch as these files can be massive."
    echo ""
	exit 0
fi

assembly=$1
hic_R1=$2
hic_R2=$3
out_prefix=$4
tmpdir=$5

PREPROC_BAM="${out_prefix}_preproc.bam"
NODUPS_BAM="${out_prefix}_nodups.bam"
NODUPS_PAIRS="${out_prefix}_nodups.pairs"

if [[ ${assembly} == *.fa ]]; then
	asm_prefix=`echo "${assembly}" | sed 's/.fa//g'`
elif [[	${assembly} == *.fasta ]]; then
	asm_prefix=`echo "${assembly}" | sed 's/.fasta//g'`
else
	echo "-> Assembly not in recognizable format (i.e., .fa or .fasta). Exiting."
fi

if [ ! -f "${assembly}.bwt" ]; then
        echo "-> BWA index not found. Indexing..."
        bwa index ${assembly}
        echo "-> Done."
else
    	echo "-> BWA index found."
fi

echo "-> Beginning alignment of Hi-C reads..."
bwa mem -SP5M -t 36 "${assembly}" "${hic_R1}" "${hic_R2}" | samtools view -bh -o "${PREPROC_BAM}"


if [ ! -f "${assembly}.fai" ]; then
	echo "-> Faidx not found. Indexing..."
	samtools faidx ${assembly}
	echo "-> Done."
	echo "-> Creating chrom.sizes file..."
	cut -f1-2 "${assembly}.fai" > "${asm_prefix}.chrom.sizes"
	echo "-> Done."
else
	echo "-> Faidx found."
	if [ ! -f "${asm_prefix}.chrom.sizes" ]; then
		echo "-> Chrom.sizes file not found. Creating..."
		cut -f1-2 "${assembly}.fai" > "${asm_prefix}.chrom.sizes"
		echo "-> Done."
	else
		echo "Chrom.sizes file found."
	fi
fi

CHROM_SIZES="${asm_prefix}.chrom.sizes"

echo "-> Beginning Pairtools pipeline."

pairtools parse --chroms-path "${CHROM_SIZES}" "${PREPROC_BAM}"
} | {
pairtools sort --nproc 24 --memory 350G --tmpdir ${tmpdir}
} | {
pairtools dedup
} | {
pairtools split --output-pairs ${NODUPS_PAIRS} --output-sam ${NODUPS_BAM}
}

if [ $? -eq 0 ]; then
	rm ${PREPROC_BAM}
	echo "-> Done."
	exit 0
else
	echo "-> Pairtools failed. Exiting."
	exit 1
fi
