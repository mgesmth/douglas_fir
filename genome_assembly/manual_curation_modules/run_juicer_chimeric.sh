#!/bin/bash
#SBATCH -J run_juicer_chimeric
#SBATCH -p general
#SBATCH -q general
#SBATCH -c 6
#SBATCH -n 1
#SBATCH --mem=20G
#SBATCH --array=[0-299]
#SBATCH -o %x.%A.%a.out
#SBATCH -e %x.%A.%a.err

set -e

echo "[M]: Host Name: `hostname`"

module load samtools/1.19
module load bwa/0.7.17
module load java-sdk/1.8.0_92

home=/home/FCAM/msmith
core=/core/projects/EBP/smith
scratch=/scratch/msmith
gid="intdf137"
site="Arima"
threads=6
jd=${core}/juicer_formanualcur

export SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID

echo "`date`:[M]: Beginning juicer run."
cd ${jd}
#Okay - now run juicer (CPU version, modified for my SLURM
#The SLURM script distributed by Aiden Lab did not work for my server
#This script is cut off after chimeric handling, and also has some additional parsing in the
#arugments section to allow the behaviour I wanted.

#i have a copy of this modified juicer script in this repo.

${jd}/scripts/juicer_justchimeric.sh -f --assembly -g "$gid" -d "${jd}/work/intdf137" -s "$site" -S chimeric \
-p references/intdf137.chrom.sizes -y restriction_sites/intdf137_Arima.txt \
-z references/interior_primary_final.fa -D "$jd" -t "$threads"

echo "`date`:[M]: Juicer chimeric handling complete! Bye."

