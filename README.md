# DLBCL Tumor Evolution

This repository contains all code used in the statistical analysis of sequencing data from serial biopsies of patients with relapsed/refractory diffuse large B-cell lymphoma. All R packages used in the analysis can be installed reproducibly with `renv`, using any R version 4.1.*. 

## Sequencing data analysis

Sequencing files were processed outside of this repository using modules from [LCR-Modules](https://github.com/LCR-BCCRC/lcr-modules) as follows: 

| Sequencing Type        | Output                             | Module(s)                |
|------------------------|------------------------------------|--------------------------|
| genome, exome, LySeqST | Somatic variant calls              | slms_3-1.0, vcf2maf-1.3  |
| genome                 | Somatic structural variants        | svar_master-1.0          |
| genome                 | Copy number                        | battenberg-1.1           |
| exome                  | Copy number                        | strelka-1.0              |
| genome                 | Subclonal structure and phylogeny  | pyclone-1.0              |
| genome, exome, LySeqST | Quality control (coverage metrics) | qc-1.0                   |
| genome, exome, LySeqST | LymphGen classifications           | lymphgen-1.0             |
| RNAseq                 | Immunoglobulin rearrangements      | bam2fastq-1.2, mixcr-1.2 |


