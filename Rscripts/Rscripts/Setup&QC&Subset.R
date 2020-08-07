#Seurat R script

#First we determine the required packages

libraries <- c("Seurat", "dplyr", "patchwork")

#Load the required packages
sapply(libraries, require, character.only=TRUE)

#Get the working directory
#Here we should pass the workign directory variable from PERL
#The software should ask for where they want to store the data

seurat.object <- setup_qc <- function(input_directory, output_directory, project, min_cells, min_features, 
                     min_nFeature_subset, max_nFeature_subset, max_percent_mt, top_features){
  
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
  
  pdf("violin_plot_mitocount_pbmc1k.pdf"); #This generates the file
  print(VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3));
  dev.off(); #This closes the graph
  
  #Also, we create a scatter plot of the features
  
  scatter_plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt");
  scatter_plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 ="nFeature_RNA");
  
  #After we generate the plots we represent them separately and together
  
  pdf("scatter_plot_features_pbmc1k.pdf");
  print(scatter_plot1)
  print(scatter_plot2)
  print(scatter_plot1 + scatter_plot2)
  dev.off();
  
  #And now we determine the working subset with the features selected
  #by the user. Here we will use the ones by defect, which are the following
  #But we have to substitute the values for those determined by the user
  
  pbmc <- subset(pbmc, subset = nFeature_RNA > min_nFeature_subset & nFeature_RNA < max_nFeature_subset & percent.mt < max_percent_mt);
  
  #Now we normalize the data from the pbmc subset using LogNormalize method with standard
  #Scale factor of 1000
  
  print ("Initializing data normalization...")
  pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  
  #We obtain the top 2000 variable features and select a subset of them defined by 
  #the user with the varaible top_features
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nFeatures = 2000)
  
  topFeatures <- head(VariableFeatures(pbmc), top_features)
  
  #Now we generate the plots for each of the data and we print them to a 
  #PDF file
  
  VFplot <- VariableFeaturePlot(pbmc)
  LabelPointPlot <- LabelPoints(plot = VFplot, points = topFeatures, repel=TRUE)
  
  pdf("VariableFeaturePlot.pdf")
  print(VFplot)
  print(LabelPointPlot)
  dev.off()
  
  #We get the names for all genes to use them as feature list for the data scaling part
  #By defect, we regress by mitochondrial dna percentage 
  
  print("Scaling data...")
  all.genes <- rownames(pbmc)
  pbmc <- ScaleData(pbmc,features = all.genes, vars.to.regress = "percent.mt")
  
  #And now we perform the dimensionality reduction by PCA
  
  print("Running dimensional reduction")
  pbmc <- RunPCA(pbmc, features=VariableFeatures(object=pbmc))
  
  #We print the PCA results to a TXT file
  
  print(pbmc[["pca"]],dims = 1:5, nFeatures=5)
  pca_dimensionality <- print(pbmc[["pca"]], dims = 1:5, nFeatures=5)
  #write.table(pbmc[["pca"]], file = "PCA_data_Seurat.txt", sep="")
  #write.table(pca_dimensionality, file = "PCA_reduced_dims.txt", sep="")
  
  #And we print the output graphs
  
  pdf("VizDimLoadings.pdf")
  print(VizDimLoadings(pbmc, dims = 1:2, reduction="pca"))
  dev.off()
  pdf("DimPlot.pdf")
  print(DimPlot(pbmc, reduction="pca"))
  dev.off()
  pdf("DimheatMap.pdf")
  print(DimHeatmap(pbmc, dims=1, cells=500, balanced=TRUE))
  print(DimHeatmap(pbmc, dims=1:15, cells=500, balanced = TRUE))
  dev.off()
  
  #And then we determine the dimensionality of the data set with JackStraw
  pbmc <- JackStraw(pbmc, num.replicate=100)
  pbmc <- ScoreJackStraw(pbmc, dims = 1:20) #We need a variable here to tell the amount of dimensions
  
  #And print the resulting plot
  pdf("JackStrawPlot.pdf")
  print(JackStrawPlot(pbmc, dims = 1:15))
  dev.off()
  
  pdf("ElbowPlot.pdf")
  print(ElbowPlot(pbmc))
  dev.off()
  
  #### HERE WE NEED TO ADD AN ALGORYTHM TO SELECT AUTOMATICALLY
  #### THE BEST JACK STRAW DIMENSIONS BUT WE ARE GOING TO KEEP
  #### THIS ORIGINAL ONE FOR THE TESTING OF THE SCRIPT
  
  ## CELL CLUSTERING
  
  #First we use the option find neighbors to find the different
  #Clusters on the chosen dimensions
  
  pbmc <- FindNeighbors(pbmc, dims = 1:10)
  pbmc <- FindClusters(pbmc, resolution = 0.5)
  
  ## NON LINEAR DIMENSIONAL REDUCTION
  
  #Now we run the dimensional reduction via umap
  
  pbmc <- RunUMAP(pbmc, dims = 1:10)
  
  #And we plot this
  
  pdf("Umap_plot.pdf")
  print(DimPlot(pbmc, reduction="umap"))
  dev.off()
  
  #And here we save the object 
  #We can use it to load it into other projects
  
  saveRDS(pbmc, file = "pbmc_tutorial.rds")
  
  ## CLUSTERING FINDING
  
  #Now we find the markers for every cluster compared to all remaining
  #cells, report only the positive ones
  
  pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
  
  
  return(pbmc) #We return the object to the parental function
}
# THIS SCRIPT ENDS HERE 
