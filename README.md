## Download Data
Commands recopilated in [download.sh]() (CAMBIAR PATH A DOWNLOAD.SH)
### 1. Download SRA files
```
prefetch SRRxxxxxxxx & #number obtained from GEO
```
### 2. Download fastq files
```
fastq-dump --split-files SRRxxxxxxx &
```

### 3. Obtain matrixes from each sample
Previously, we have to change the name to the correct format: 
```
mv SRRxxxxxxx_1.fastq <sample>_S1_L001_R1_001.fastq
mv SRRxxxxxxx_2.fastq <sample>_S1_L001_R2_001.fastq
```
or
```
mv SRRxxxxxxx_1.fastq <sample>_S1_L001_I1_001.fastq
mv SRRxxxxxxx_2.fastq <sample>_S1_L001_R1_001.fastq
mv SRRxxxxxxx_3.fastq <sample>_S1_L001_R2_001.fastq
```
Then: 
```
#If the data were obtained with 10x Genomics: 
cellranger count --id=<sample> --transcriptome=/CEPH/shared/databases/cellranger/refdata-gex-GRCh38-2020-A --fastqs=<fastqs directory> --include-introns=true --sample=<sample> --localcores=30 --localmem=64 &
```

## Preparation of data
Scripts recopilated in [execution.sh](https://github.com/mariavam/ADSC/blob/e813a092c0852a398b6b3126c708beb337571a48/execution.sh)

For their execution, all of them reference the follow documents:
  + [libraries.sh]() contain all libraries used
  + [parameters.sh]() contain the parameters
