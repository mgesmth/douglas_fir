#!/bin/bash


module load fastqc/0.12.1
module load MultiQC/1.10.1 

cd ${outdir}


fastqc -o ./ -t 12 -d ${read1} ${reads2}
multiqc .
