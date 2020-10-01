## SCRIPT: ANALYSIS OF scRNAseq Data

libraries <- c("Seurat", "dplyr", "patchwork")
sapply(libraries, require, character.only=TRUE)

setup_qc <- function(input_directory, output_directory, project_name, min_cells, min_features, 
                     min_nFeature_subset, max_nFeature_subset, max_percent_mt){
  
  setwd(output_directory)
  
  #Load the PBMC dataset with the data from Cellranger
  #We should add the directory to the script from perl
  project.data <- Read10X(data.dir = input_directory);
  
  #Initialize the Raw Seurat Object
  
  print ("Creating Seurat Object...")
  project_name <- CreateSeuratObject(counts = project.data, project = project_name, min.cells = min_cells, min.features = min_features);
  
  #After the object has been created we initiate the standard pre-processing
  #in which we obtain the mitochondrial DNA
  
  project_name[["percent.mt"]] <- PercentageFeatureSet(project_name, pattern = "^MT-");
  
  #Now we export the graphs with different visualization
  
  #We create a violin plot of the data and export it to pdf
  jpeg(filename = "Violin_plot.jpg", units ="px", width = 1024, height = 1024)  #Opens the picture
  print(VlnPlot(project_name, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3));
  dev.off(); #This closes the graph
  
  #Also, we create a scatter plot of the features
  
  scatter_plot1 <- FeatureScatter(project_name, feature1 = "nCount_RNA", feature2 = "percent.mt");
  scatter_plot2 <- FeatureScatter(project_name, feature1 = "nCount_RNA", feature2 ="nFeature_RNA");
  
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
  
  project_name <- subset(project_name, subset = nFeature_RNA > min_nFeature_subset & nFeature_RNA < max_nFeature_subset & percent.mt < max_percent_mt);
  
  #And we return the object to the function
  
  return(project_name); 
}

#The second function we are going to design is the one implicated on
#data normalization and feature selection 

print("Setup and Quality control function... Loaded!")

Normalization <- function(parse_object, norm_method, scfactor,
                         selection_method, number_features, top_features){
  
  function_object <- parse_object
  
  #We firstly have to normalize the data
  
  function_object <- NormalizeData(function_object, normalization.method = norm_method,
                        scale.factor = scfactor)
  
  #After data normalization we have to identify the highly variable features
  
  function_object <- FindVariableFeatures(function_object, selection.method = selection_method,
                               nfeatures = number_features)
  
  #We identify the top N variable features, which are the ones indicated by the
  #user in the function
  
  topNfeatures <- head(VariableFeatures(function_object), top_features)
  
  #And we print the top N features we have gathered
  
  NFeatures_plot1 <- VariableFeaturePlot(function_object)
  NFeatures_plot2 <- LabelPoints(plot = NFeatures_plot1, points = topNfeatures, repel = TRUE)
  jpeg(filename = "Variable_featurePlot.jpg", units ="px", width = 2048, height = 1024)  #Opens the picture
  print(NFeatures_plot1 + NFeatures_plot2)
  dev.off()
  
  
  return(function_object)
  
}

print("Normalization function... loaded!")

## Third function will be in charge of PCA and determination of the 
## dimensionality of the dataset

PCA_and_dimensionality <- function(parse_object, dimensions_on_text,
                                   features_on_text, VDL_dims, DHM_dims,
                                   DHM_cells){
  
  function_object <- parse_object
  
  #We scale the data 
  
  all.genes <- rownames(function_object)
  function_object <- ScaleData(function_object, features = all.genes)
  
  #The first thing is to run the PCA reduction, with the features as the 
  #variable features from the original object
  
  function_object <- RunPCA(function_object, features = VariableFeatures(
    object = function_object))
  
  #After doing this reduction, we will visualize the PCA in different ways
  #The first is by text, where we will print to a .txt file the output of the
  #desired dimensions the user tells with the variable dimensions_on_text
  
  capture.output(print(function_object[["pca"]], dims = 1:dimensions_on_text, 
                       nfeatures = features_on_text), file = "PCA_dims.txt")
  
  #Now we get different plots and print them to PDF files (they're vectorial)
  #To print them we will use a if-else loop in which each plot will be plotted
  #In a different page of the PDF file. This way we will have a full-size plot
  #That the user will use accordingly. This will be applied for both 
  #VizDimLoadings and DimHeatMap plots where we can have several ones. The one
  #for DimPlot doesn't need a loop as its a single graph
  
  pdf(file = "VizDimLoadingPlot.pdf")
  if(VDL_dims > 1){
    i = 1
    while(i <= VDL_dims){
      #plot_name <- paste("VizDimLoading_dimension_", i, ".jpeg", sep = "")
      #jpeg(file = plot_name, width = 1024, height = 1024)
      plotA <- VizDimLoadings(function_object, dims = i, reduction = "pca")
      print(plotA)
      #dev.off()
      i <- i+1
    }
  }
  else if(VDL_dims == 1){
    plotA <- VizDimLoadings(function_object, dims = VDL_dims, reduction = "pca")
    print(plotA)
  }
  else{
    print("The value for dimensions on VDL_dims is not valid, exiting...")
    break
  }
  dev.off()
  
  pdf(file = "DimPlot.pdf")
  graph2 <- DimPlot(function_object, reduction = "pca")
  print(graph2)
  dev.off()
  
  pdf(file = "DimHeatMaps.pdf")
  plot_DHM <- DimHeatmap(function_object, dims = 1:DHM_dims, cells = DHM_cells, balanced = TRUE)
  print(plot_DHM)
  dev.off()
  
  return(function_object)
}

print("PCA calculation function... loaded!")

dimensionality <- function(parse_object, JS_replicates, JSscore_dims,
                           JSPlot_dims){

  #First we assign the object we get to a variable in the subroutine
  
  function_object <- parse_object

  #After this, we can operate with the seurat object containing all the previous data

  function_object <- JackStraw(function_object, num.replicate = JS_replicates)
  
  #We perform a JackStraw analysis on the object with a number of replicates stablished
  #by the user, with a default value of 100 replicates. After that, we analyze the score 
  #of the Jack Straw dimension analysis. We use, by default, the dimensions until 20, but
  #the user can modify this value accordingly
  
  function_object <- ScoreJackStraw(function_object, dims  = 1:JSscore_dims)
  
  #Now we plot different plots
  
  JSPlot <- JackStrawPlot(function_object, dims = 1:JSPlot_dims)
  ElPlot <- ElbowPlot(function_object)
  
  pdf(file = "JackStrawPlot.pdf")
  print(JSPlot)
  dev.off()
  
  pdf(file = "ElbowPlot.pdf")
  print(ElPlot)
  dev.off()
  
  #After doing the plots, we will automatically select the amount of dimensions
  #that are representative of the PC variation. For this we will use the values of 
  #stdev obtained in the Elbow Plot, calculate the median value for them and automatically
  #select all those dimensions over this. 
  
  return(function_object)
}

print("Dimensionality calculation function... loaded!")

clustering_and_NLR <- function(parse_object, PC_chosen, cluster_resolution,
                      reduction_method){
  
  function_object <- parse_object
  
  #We first find the neighbors of the clusters. The dimensions are the ones 
  #described by the user previously based on the information obtained by JackStraw

  function_object <- FindNeighbors(function_object, dims = 1:PC_chosen)
  
  #And after this we obtain the clusters. The value needed is the resolution,
  #that will be provided with the variable cluster_resolution
  
  function_object <- FindClusters(function_object, resolution = cluster_resolution)
  
  #With this information we can see the IDs of the first N cells asked by the user
  #in the cluster
  
  options(max.print= 99999999) #With this we avoid the limit for printing
  
  capture.output(Idents(function_object), file = "Cluster_IDs.txt")
  
  #After doing this clustering we perform the non-linear reduction of the data
  #using the same dimensions we already chose
  
  function_object <- RunUMAP(function_object, dims = 1:PC_chosen)
  
  #And we generate the DimPlot
  
  DPlot <- DimPlot(function_object, reduction = reduction_method)
  
  pdf(file = "DimPlot.pdf")
  print(DPlot)
  dev.off()
  
  return(function_object)
}

print("Clustering function... loaded!")

#After doing all this tasks we will have several .Rds files with the different steps
#and the date they were performed. We encourage the users to open any of the files and
#explore all the possiblities there


