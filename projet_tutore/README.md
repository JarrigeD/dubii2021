## Introduction
#### This repository contains my personal training project for DUBII. I worked on it under the guidance of [clairetn](https://github.com/clairetn) at [I2BC](https://www.i2bc.paris-saclay.fr/spip.php?article391).

The subject was: 
# Study of nuclear and organellar genomes interactions in the green microalga _Chlamydomonas reinhardtii_.

It was based on public data from :

**[Strenkert D, Schmollinger S, Gallaher SD, et al. Multiomics resolution of molecular events during a day in the life of Chlamydomonas. Proc Natl Acad Sci U S A. 2019;116(6):2374-2383. doi:10.1073/pnas.1815238116](https://www.pnas.org/content/116/6/2374)**

Where synchronised cultures of _Chlamydomonas reinhardtii_ (strain CC-5390) were grown in photobioreactors on a 12h light / 12h dark cycle.
Samples were collected at regular intervals and subjected to a transcriptomics analysis by RNA-Seq and proteomics analysis by LC-MS/MS.

The RNA-seq data was collected from: https://www.ebi.ac.uk/ena/browser/view/PRJNA445880

The protein abundance table from: https://www.pnas.org/content/suppl/2019/01/17/1815238116.DCSupplemental

## Contents

- snakemake workflow: to prepare the RNA-seq data from raw fastqc to gene counts table
- R script: to analyse the transcriptomic data as well as the proteomic
- mini RNA-seq dataset to test the workflow: 5 mini-fastq.gz files containing the first 1000 reads from the original RNA libraries ([Strenkert D, Schmollinger S, Gallaher SD, et al. Multiomics resolution of molecular events during a day in the life of Chlamydomonas. Proc Natl Acad Sci U S A. 2019;116(6):2374-2383. doi:10.1073/pnas.1815238116](https://www.pnas.org/content/116/6/2374))
- metadata file: with information on the miniRNA-seq dataset

## How to use the workflow and scripts:
(https://snakemake.github.io/) [![Snakemake](https://img.shields.io/badge/snakemake-≥5.6.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
### Required:
- [snakemake]
- slurm-drmaa
- STAR
- featureCounts
- R
- IN PROGRESS
- 
You can use the mini dataset to test the snakemake workflow, change the parameters in the config file before running. Build the project architecture with the various directories first. 

### To run workflow on slurm cluster:

`module load snakemake`

`module load slurm-drmaa`

`snakemake --drmaa  --jobs={24} --snakefile {scripts/rna_workflow.smk} --configfile {scripts/config.yml}`


IN PROGRESS
