# FIRST SCRIPT FOR SEURAT ANALYSIS

#The first thing we have to do is to load all the required dependencies

libraries <- c("Seurat", "dplyr", "patchwork")
sapply(libraries, require, character.only=TRUE)


setup_qc <- function(input_directory, output_directory, project, min_cells, min_features, 
                                      min_nFeature_subset, max_nFeature_subset, max_percent_mt){
  
  setwd(output_directory)
  
  #Load the PBMC dataset with the data from Cellranger
  #We should add the directory to the script from perl
  pbmc.data <- Read10X(data.dir = input_directory);
  
  #Initialize the Raw Seurat Object
  
  print ("Creating Seurat Object...")
  pbmc <- CreateSeuratObject(counts = pbmc.data, project = project, min.cells = min_cells, min.features = min_features);
  
  #After the object has been created we initiate the standard pre-processing
  #in which we obtain the mitochondrial DNA
  
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-");
  
  #Now we export the graphs with different visualization
  
  #We create a violin plot of the data and export it to pdf
  jpeg(filename = "Violin_plot.jpg", units ="px", width = 1024, height = 1024)  #Opens the picture
  print(VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3));
  dev.off(); #This closes the graph
  
  #Also, we create a scatter plot of the features
  
  scatter_plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt");
  scatter_plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 ="nFeature_RNA");
  
  #After we generate the plots we represent them separately and together
  
  jpeg(filename = "ScatterPlot_nCountRNA_percentmt.jpg", units ="px", width = 1024, height = 1024)  #Opens the picture
  print(scatter_plot1)
  dev.off()
  jpeg(filename = "ScatterPlot_nCountRNA_nFeatureRNA.jpg", units ="px", width = 1024, height = 1024)  #Opens the picture
  print(scatter_plot2)
  dev.off()
  jpeg(filename = "ScatterPlot_together.jpg", units ="px", width = 2048, height = 1024)  #Opens the picture
  print(scatter_plot1 + scatter_plot2)
  dev.off();
  
  #And now we determine the working subset with the features selected
  #by the user. Here we will use the ones by defect, which are the following
  #But we have to substitute the values for those determined by the user
  
  pbmc <- subset(pbmc, subset = nFeature_RNA > min_nFeature_subset & nFeature_RNA < max_nFeature_subset & percent.mt < max_percent_mt);
  
  #And we return the object to the function
  
  return(pbmc); 
}

#We gather the variables from the script

project_name <- "pbmc3k"
input <- "C:/Users/gvast/OneDrive/Documentos/TFM/dataset_pbmc3k/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19" #esto vendría del script PERL
output <- "C:/Users/gvast/OneDrive/Escritorio/Prueba3"
minimum_cells <- 3
minimum_features <- 200
minimum_subset_features <- 200
maximum_subset_features <- 2500
maximum_subset_mt <- 5

#Now we generate the object, using the structure "SeuratObject_date"

date <- Sys.Date()
name <- paste("SeuratObject", project_name, date, sep="_")

#After executing the scritp we use the function assign to create an object
#With the name stored in the variable name and the result of the object 
#generated in the function

assign(name, setup_qc(input, output, project_name, minimum_cells,
                      minimum_features, minimum_subset_features,
                      maximum_subset_features, maximum_subset_mt))

#And after this we would save the object into the workign directory as
#the name specified with the date of the analysis

saveRDS(toString(name), file=paste(toString(name),"QCScript.rds",sep=""))

     