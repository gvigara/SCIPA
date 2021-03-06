# SCIPA: Single Cell Integrated-Pipeline for Analysis

## Introduction

Recent advance in resolution and capacity of single cell RNA seq has allowed for an increasing use of this technique. The information avaiable on a single run
is massive, allowing from the identification of multiple cell populations on tissues to the identification and isolation of transcriptional changes on a group of
cells specifically targeted by a virus. This rise on the use of scRNAseq, alongside with the improvement in the used technology, has lead to an increasing amount of
experiments as well as an increase in the data contained in those experiments. As the amount of cells and complex of the experiment increases, so does the weight of the
analysis files, and with it the time required for their analysis. 

SCIPA stands for Single Cell Integrated-Pipeline for Analysis. This is a pipeline developed using two powerful programming languages, Perl and R, oriented to the
automation of analysis on single cell RNA seq dataset, based on 10X chromium technology (both V2 and V3). For this means, SCIPA is build around a main program that
acts as a principal pipeline communicating several Perl modules, written for different specific tasks, alongside with an R script written for automatic analysis
of the expression matrices obtained throught the process.

<p align="center">
  <img src="https://github.com/gvigara/SCIPA/tree/master/img/Image_1.jpg">
</p>

## Download

To download the software, please download the repo and uncompress it on the desired output folder:

```
$ curl -L -O https://github.com/gvigara/SCIPA
$ unzip SCIPA
$ cd SCIPA
```

When in the folder you will see several files. The config.ini file will be the configuration file where all the parameters are set. Please follow inside instructions for each parameter. 

## Installation and usage

First use the dependencies script to install all necessary dependencies. This script will ask for the super-user password in order to install them. If the user wants
to install manually the dependencies just open the file on a text editor (g.e Kate) and execute all the commands manually. 

```
$ chmod +x dependencies.sh
$ sudo ./dependencies.sh
```

Dependencies include: 

 - curl
 - fastqc
 - python3-pip
 - multiqc (from pip)
 - libcurl4-openssl-dev
 - libssl-dev
 - r-base-core (for those who doesn't have R already installed, you can also download it form [R](https://www.r-project.org/))
 - r-cran-httr
 - r-cran-seurat
 - r-cran-dplyr
 - r-cran-devtools
 - r-bioc-limma
 
In the case that there is any problem during R execution you may have to add your user to the staff group:

`$ sudo useradd -a -G staff user`

And after that give reading permisions to the library folder:

`$ sudo chmod o+w /usr/local/lib/R/site-library`

A major need is the download of cellranger form 10X genomics webpage. The user has to uncompress cellranger on the wanted directory and then locate this directory on the config.ini file. Cellranger can be downloaded from this [link](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest)
The output of each analysis will be generated on a dedicated folder on the working directory designated on the config.ini file. The Working directory is defined as the directory where the FASTQ files are located, and will be condisered the root directory for the job. 

## Whole analysis process

This is meant to be an example on how to analyze a small sample, like PBMC 1k V3, downloaded from [10X genomics webpage](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_v3). This .zip file has to be uncompressed on a directory that will be the root directory for our work (designed as Working directory on the config.ini file). 

Afterwards we have to set all the parameters on the config.ini file, stablishing the following: 

- Directory where Cellranger is located
- Type of organism (Human/Mouse). In case of Human selection, choose 1 for GRCh38 or 2 for HG19. 
- 
