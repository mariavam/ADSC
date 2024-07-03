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
Scripts recopilated in [execution.sh](https://github.com/mariavam/ADSC/blob/e813a092c0852a398b6b3126c708beb337571a48/execution.sh).

This Analysis is performed with single-cell RNA seq datasets from 5 different datasets. A Differential Expression Analysis of these dataset, and compare them with the results obtained from in-vitro obtained samples. 

For their execution, all of them reference the follow documents:
  + [libraries.r](https://github.com/mariavam/ADSC/blob/master/libraries.r) contain all libraries used
  + [parameters.r](https://github.com/mariavam/ADSC/blob/master/parameters.r) contain the parameters

Each dataset were prepared with the files [prep_<papername>.r]. There is one file per paper (Grubman, Leng, Otero, Alsema and Mathys).
Then they were filtered with the code of [filter.r](https://github.com/mariavam/ADSC/blob/master/filter.r)

## Analysis
Then, the Post Mortem analysis was performed with the code of [PM_Analysis.Rmd](https://github.com/mariavam/ADSC/blob/master/PM_Analysis.Rmd) and [PM_DEA.Rmd](https://github.com/mariavam/ADSC/blob/master/PM_DEA.Rmd), that contains the differential expression analysis. 

The comparison was performed with the code from [Comparison_IV_PM.Rmd](https://github.com/mariavam/ADSC/blob/master/Comparison_IV_PM.Rmd).

## Session info
R version 4.3.3 (into account --> R.studio version 4.2.1)

