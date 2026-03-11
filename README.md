# douglas_fir
Code for interior Douglas-fir (_Pseudotsuga menziesii_ var. _glauca_) de novo assembly project.


# Genome Assembly

## Polishing and Manual Curation

Upon curating Hi-C maps with Juicer on initial scaffolded assembly, several major misassemblies were noted. To correct these, we used Juicer and the 3D-DNA pipeline to computationally and, using JBAT tools, manually correct these misassemblies.

The Juicer program has a SLURM version. For my particular server, the script is not compatible, and I found it more effective to modify the CPU version of Juicer to suit my needs. Additionally, Juicer nativelty relies on SAM alignment files, which may be feasible for some, but given the size of the Douglas-fir genome, I needed store alignments in BAM format. So I modified the scripts to store alignmenmts as BAMs, and then pass decompressed SAM files to Juicer commands through `samtools view` and a pipe. Different versions of modified Juicer scripts are found here: `genome_assembly/manual_curation_modules/juicer_*.sh`

For the 3D-DNA pipeline, due to the size and low complexity of the Douglas-fir genome (~70% repetitive), we increased splitting, editing and polishing parameters to make the editing less sensitive and correcting false misjoins.

```

#prep Hi-C data, draft assembly and environment for Juicer
./genome_assembly/manual_curation_modules/prep_forjuicer.sh

#align split HiC files to draft assembly in parallel
./genome_assembly/manual_curation_modules/align_hic.sh

#handle chimeric alignments in parallel
./genome_assembly/manual_curation_modules/run_juicer_chimeric.sh

#run Juicer in assembly mode from merge (post-chimeric, pre-merge)
./genome_assembly/manual_curation_modules/run_juicer_merge.sh

#run 3D-DNA pipeline with two rounds of iterative misjoin correction
./genome_assembly/manual_curation_modules/3DDNA.sh

```
