## SEURAT script for Data analysis using SCTransform

#DEPENDENCY CHECK SCRIPT

#Seurat software is dependent on the following packages:
#Seurat (httr and plotly are dependecies of seurat)
#GGplot2 (installed)
#Dplyr
#Patchwork


dependencies <- c("RCurl", "httr", "plotly", "Seurat","patchwork","dplyr") #All required pacakges

#We create the function to install all the packages if missing

print("Checking dependencies...")``

use.packages <- function(pkgs){
  new.pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])] #check the list of installed packages
  if(length(new.pkgs)) #If not installed
    install.packages(new.pkgs, dependencies = TRUE)
  sapply(pkgs, require, character.only = TRUE) #This loads the packages
}

#If they are not missing, the use of the function use.packages will load them

use.packages(dependencies) #Now we use the function with the packages
