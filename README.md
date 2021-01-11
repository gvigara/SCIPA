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

![Image](https://github.com/gvigara/SCIPA/tree/master/img/SCIPA_diagram.jpg)

## Installation and usage

First use the dependencies script to install all necessary dependencies. This script will ask for the super-user password in order to install them. If the user wants
to install manually the dependencies just open the file on a text editor (g.e Kate) and execute all the commands manually. 

`$ chmod +x dependencies.sh`

`$ sudo ./dependencies.sh`

Dependencies include: 

 - curl
 - fastqc
 - python3-pip
 - multiqc (from pip)
 - libcurl4-openssl-dev
 - libssl-dev
 - r-base-core (for those who doesn't have R already installed)
 - r-cran-httr
 - r-cran-seurat
 - r-cran-dplyr
 - r-cran-devtools
 - r-bioc-limma
 
In the case that there is any problem during R execution you may have to add your user to the staff group:

`$ sudo useradd -a -G staff user`

And after that give reading permisions to the library folder:

`$ sudo chmod o+w /usr/local/lib/R/site-library`
