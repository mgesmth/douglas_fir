# Douglas-fir Genomics Pipeline

Code is being added as the project progresses.

## Assembly Quality Assessment

### Statistics 

QUAST/, meryl/, merqury/, GenomeScope/ and BUSCO/ were used to gather quality statistics for both interior and coastal Douglas-fir genome assemblies as appropriate. 

```sh
#The code used to build the intDF011 assembly (not CBP):
./intDF011_pipeline.sh

#Assembly quality statistics for CBP and coastal assemblies, as well as Hi-C contact maps for interior assemblies:
./assembly_quality_assessment.sh
```
