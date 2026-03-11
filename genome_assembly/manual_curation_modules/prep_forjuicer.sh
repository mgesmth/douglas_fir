#!/bin/bash
#SBATCH -J prep_forjuicer
#SBATCH -p himem
#SBATCH -q himem
#SBATCH -c 24
#SBATCH --mem=128G
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err

set -e

date
echo "[M]: Host Name: `host name`"

module load bwa/0.7.17
module load samtools/1.20
module load python/3.8.1
module load seqkit/2.10.0
home=/home/FCAM/msmith
core=/core/projects/EBP/smith
scratch=/scratch/msmith
juicedir=${core}/juicer_formanualcur
export PATH="${juicedir}/scripts:$PATH"
asm=${juicedir}/references/interior_primary_final.fa
asm_name=$(basename ${asm})
gid="intdf137"
enzyme="Arima"
workdir=${juicedir}/work/${gid}

if [[ ! -d ${workdir} ]] ; then
  mkdir ${workdir}
fi

R1=${home}/hiC_data/allhiC_R1.fastq.gz
R2=${home}/hiC_data/allhiC_R2.fastq.gz
hic_split=${workdir}/fastq
hic_bams=${workdir}/splits
splitN=300

if [[ ! -d "$hic_split" ]] ; then
  mkdir ${hic_split}
fi


echo -e "\n[M]: Splitting Hi-C data\n"

seqkit split2 -f -1 "$R1" -2 "$R2" -p "$splitN" -O "$hic_split" 
cd $hic_split
ls allhiC_R1*.fastq.gz > fastqs.txt

date=$(date)
echo -e "\n${date}:[M]: HiC data split. Moving onto site_positions file."

cd ${juicedir}/restriction_sites
python ${juicedir}/scripts/generate_site_positions.py "$enzyme" "$gid" "$asm"

date=$(date)
echo -e "\n${date}:[M]: Site positions files generated."

#check if there's bwa index has been run on the assembly already; if not, run it
if [[ ! -f ${asm}.bwt ]] ; then
  echo -e "${date}:[M]: No BWA index detected. Running bwa index."
  bwa index ${asm}
  date=$(date)
  echo -e "\n${date}:[M]: BWA index created."
else
  echo -e "${date}:[M]: BWA index found."
fi

echo -e "\n${date}:[M]: Complete. Bye."
