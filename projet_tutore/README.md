## Introduction
### This repository contains my personal training project for DUBII. I worked on it under the guidance of [clairetn](https://github.com/clairetn) at [I2BC](https://www.i2bc.paris-saclay.fr/spip.php?article391).

The subject was: 
# Study of nuclear and organellar genomes interactions in the green microalga _Chlamydomonas reinhardtii_.

It was based on public data from :

**[Strenkert D, Schmollinger S, Gallaher SD, et al. Multiomics resolution of molecular events during a day in the life of Chlamydomonas. Proc Natl Acad Sci U S A. 2019;116(6):2374-2383. doi:10.1073/pnas.1815238116](https://www.pnas.org/content/116/6/2374)**

Where synchronised cultures of _Chlamydomonas reinhardtii_ (strain CC-5390) were grown in photobioreactors on a 12h light / 12h dark cycle.
Samples were collected at regular intervals and subjected to a transcriptomics analysis by RNA-Seq and proteomics analysis by LC-MS/MS.

The RNA-seq data was collected from: https://www.ebi.ac.uk/ena/browser/view/PRJNA445880

The protein abundance table from: https://www.pnas.org/content/suppl/2019/01/17/1815238116.DCSupplemental

## Contents
You will find here:
- the main snakemake workflow I wrote to prepare the RNA-seq data from raw fastqc to gene counts table
- the R script I wrote to analyse the transcriptomic data as well as the proteomic
- a mock dataset to test the workflow

## How to use the workflow and scripts:
### Required:
- (snakemake)[https://snakemake.github.io/]
- STAR
- featureCounts
- R
- IN PROGRESS
- 
You can use the mock dataset to test the snakemake workflow 

IN PROGRESS
