#!/bin/bash
# This script will install all the dependencies
# required for the pipeline to function properly
#
# To run this file, you need to be the superuser
# so, upon download, just execute:
# chmod +x dependencies.sh
# sudo ./dependencies.sh

#First we check if the script is being executed as superuser
if(($(id -u) != 0)) #If sudo, id -u will return 0
then
    echo "This script has to be executed as superuser. Please, execute sudo ./dependecies.sh"
    exit
fi

#If the script is executed as superuser we can now begin the installation
#First we need to update all the repositories

apt-get update

#Now we will install curl (in case it's not installed)

printf "\n\n Installing curl..."
apt-get -y install curl

#And FASTQC for quality control of the readings

printf "\n\n Installing FASTQC..."
apt-get -y install fastqc

#We need pip for installing multiqc

printf "\n\n Installing python pip"
apt-get -y install python3-pip

#And now we can install multiqc from the repositories

printf "\n\n Installing multiqc..."
pip install multiqc

##  We need R installed along with several libraries for analysis

printf "\n\n Installing dependencies for R..."

apt-get -y install libcurl4-openssl-dev libssl-dev r-base-core

# We need to install Statistics::R for the Seurat analysisi

printf "\n\n Installing Statistics module for Perl..."

perl -MCPAN -e 'install Statistics::R'

# And some libraries for R

printf "\n\n Installing R libraries..."

apt-get -y install r-cran-httr r-cran-seurat r-cran-dplyr r-cran-devtools

