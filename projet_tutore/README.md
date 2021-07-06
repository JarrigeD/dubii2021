## Introduction
#### This repository contains my personal training project for DUBII. I worked on it under the guidance of [clairetn](https://github.com/clairetn) at [I2BC](https://www.i2bc.paris-saclay.fr/spip.php?article391).

The subject was: 
# Study of nuclear and organellar genomes interactions in the green microalga _Chlamydomonas reinhardtii_.

It was based on public data from :

**[Strenkert D, Schmollinger S, Gallaher SD, et al. Multiomics resolution of molecular events during a day in the life of Chlamydomonas. Proc Natl Acad Sci U S A. 2019;116(6):2374-2383. doi:10.1073/pnas.1815238116](https://www.pnas.org/content/116/6/2374)**

Where synchronised cultures of _Chlamydomonas reinhardtii_ (strain CC-5390) were grown in photobioreactors on a 12h light / 12h dark cycle.
Samples were collected at regular intervals and subjected to a transcriptomics analysis by RNA-Seq and proteomics analysis by LC-MS/MS.

The RNA-seq data was collected from: 

https://www.ebi.ac.uk/ena/browser/view/PRJNA445880

The protein abundance table from: 

https://www.pnas.org/content/suppl/2019/01/17/1815238116.DCSupplemental

## Contents

- snakemake workflow: to prepare the RNA-seq data from raw fastqc to gene counts table 
- associated R script: to analyse the transcriptomic data and extract rhythmic and differentially expressed genes
- conda environment files
- template configuration file (to edit before running the workflow)
- mini RNA-seq dataset to test the workflow: 5 mini-fastq.gz files containing the first 1000 reads from the original RNA libraries ([Strenkert D, Schmollinger S, Gallaher SD, et al. Multiomics resolution of molecular events during a day in the life of Chlamydomonas. Proc Natl Acad Sci U S A. 2019;116(6):2374-2383. doi:10.1073/pnas.1815238116](https://www.pnas.org/content/116/6/2374))
- metadata file: with information on the miniRNA-seq dataset
- html report showing part of my statistical analyses in R


## How to use the workflow and scripts:
[![Snakemake](https://img.shields.io/badge/snakemake-≥5.6.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
### Required:
- [snakemake](https://snakemake.github.io/) 
- [conda](https://docs.conda.io/projects/conda/en/latest/index.html)

You can use the mini dataset to test the snakemake workflow, change the parameters in the config file before running. Build the project architecture with the various directories first. The test data is too small to extract significant statistical results with the R script part of the workflow.

Instead, this **[html report](https://jarriged.github.io/dubii2021/projet_tutore/rhythmic_analyses)** shows a glimpse of what can be found using the statistical analysis with the data from _Strenkert et al, 2019_.

### To run workflow:

`snakemake --cores={your_pick} --snakefile {rna_workflow.smk} --configfile {config.yml} --use-conda`


IN PROGRESS
