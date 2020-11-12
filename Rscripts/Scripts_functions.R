## SCRIPT: ANALYSIS OF scRNAseq Data

libraries <- c("Seurat", "dplyr", "patchwork", "future")
sapply(libraries, require, character.only=TRUE)

setup_qc <- function(input_directory, output_directory, project_name, min_cells, min_features, 
                     min_nFeature_subset, max_nFeature_subset, max_percent_mt){
  
  #Print the variables for execution
  print("Initiating setup and QC analysis. Function parameters:")
  message("Input directory: ", input_directory)
  message("Output directory: ", output_directory)
  message("Project name: ", project_name)
  message("Minimum cells for analysis: ", min_cells)
  message("Minimum features for analysis: ", min_features)
  message("Minimum subset features: ", min_nFeature_subset)
  message("Maximum subset features: ", max_nFeature_subset)
  message("Maximum subset mitoDNA: ", max_percent_mt)               

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
  
  #If the user doesn't add any information about the values for the subset
  #Now we determine the min, max features for the subset and the max mitoDNA
  
  if(min_nFeature_subset == "." & max_nFeature_subset == "." & max_percent_mt == "."){
  minimum_subset_features = min(project_name[["nFeature_RNA"]]$nFeature_RNA)
  #As maximum subset we will use the 99% of the values, gathered under percentile 99.9
  maximum_subset_features = quantile(project_name[["nFeature_RNA"]]$nFeature_RNA, 0.999)
  #As mitochondrial value we will use the values under the percentile 95, gathering most of the values there
  maximum_mitochondral_dna_percentage = quantile(project_name[["percent.mt"]]$percent.mt, 0.95)
  print("The user didn't provide any parameters for subset generation, using the calculated ones")
  print(paste("Minimum subset features:", minimum_subset_features, sep = " "))
  print(paste("Maximum_subset_features:", maximum_subset_features, sep = " "))
  print(paste("Maximum_mitochondrial_dna_percentage:", maximum_mitochondral_dna_percentage, sep = " "))
  project_name <- subset(project_name, subset = nFeature_RNA > minimum_subset_features & nFeature_RNA < maximum_subset_features & percent.mt < maximum_mitochondral_dna_percentage)
  }
  else{
    #If the user adds information about the subset we use its values
    project_name <- subset(project_name, subset = nFeature_RNA > min_nFeature_subset & nFeature_RNA < max_nFeature_subset & percent.mt < max_percent_mt);
    
  }

    #And we return the object to the function
  
  return(project_name); 
}

#The second function we are going to design is the one implicated on
#data normalization and feature selection 

print("Setup and Quality control function... Loaded!")

Normalization <- function(parse_object, norm_method, scfactor,
                         selection_method, number_features, top_features, cores){
  
  function_object <- parse_object

  print("Initiating normalization function. Parameters:")
  message("Normalization method: ", norm_method)
  message("Scaling factor: ", scfactor)
  message("Selection method: ", selection_method)
  message("Number of features: ", number_features)
  message("Top variable features", top_features)
  
  # we set the paralelization value
  plan(multicore, workers = cores) ########### IMPORTANTE
  
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
                                   DHM_cells, cores){
  
  function_object <- parse_object
  plan(multicore, workers = cores)
  
  print("Starting PCA and dimensionality adjustment:")
  message("Dimensions on text: ", dimensions_on_text)
  message("Features on text: ", features_on_text)
  message("VizDimLoadings dimensions: ", VDL_dims)
  message("DimHeatMap dimensions: ", DHM_dims)
  message("DHM cells for analysis: ", DHM_cells)
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
                           JSPlot_dims, cores){

  #First we assign the object we get to a variable in the subroutine
  
  function_object <- parse_object
  
  

  print("Calculating dimensionality:")
  message("JackStraw replicates: ", JS_replicates)
  message("JackStrawScore dimensions: ", JSscore_dims)
  message("JackStrawPlot dimensions: ", JSPlot_dims)

  #After this, we can operate with the seurat object containing all the previous data
  
  plan(multicore, workers = cores)
  function_object <- JackStraw(function_object, num.replicate = JS_replicates, dims = JSscore_dims) #50 is the maximum amount of dimensions to be performed
  #we change all the values to 50 so we can calculate all the scores and select from them 
  
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

clustering_and_NLR <- function(parse_object, PC_pvalue, cluster_resolution,
                               reduction_method, organism, selected_tissue, cores){
  
  function_object <- parse_object
  
  plan(multicore, workers = cores)

  print("Clustering and Neighbour calculation: ")
  message("Principal componentes p-value threshold: ", PC_pvalue)
  message("Cluster resolution analysis: ", cluster_resolution)
  message("Reduction method for analysis: ", reduction_method)
  
  #First we need to choose the principal components using the algorithym:
  
  total_dimensions = length(function_object[["pca"]]@jackstraw$overall.p.values[,1]) #This way we get the total amount of PC analyzed
  i <- 1
  p <- 1 #Where we start to add elements to the list
  selected_PC <- vector() #We generete the empty vector that will store the PC chosen
  while(i <= total_dimensions){
    if(function_object[["pca"]]@jackstraw$overall.p.values[i,2] <= PC_pvalue){
      print(paste("Dimension", i, "under p-value threshold, adding to chosen PC list...", sep= " "))
      selected_PC[p] <- i
      p <- p + 1
      i <- i+1
    }
    else{
      print(paste("Dimension", i, "over p-value threshold, skipping...", sep=" "))
      i <- i+1
    }
  }
  
  capture.output(
  paste(print("Selected_PCs_are:"),
  print(selected_PC), sep = ""), file = "Selected_PCs.txt")

  
  #We first find the neighbors of the clusters. The dimensions are the ones 
  #described by the user previously based on the information obtained by JackStraw
  
  function_object <- FindNeighbors(function_object, dims = selected_PC)
  
  #And after this we obtain the clusters. The value needed is the resolution,
  #that will be provided with the variable cluster_resolution
  
  function_object <- FindClusters(function_object, resolution = cluster_resolution)
  
  #With this information we can see the IDs of the first N cells asked by the user
  #in the cluster
  
  options(max.print= 99999999) #With this we avoid the limit for printing
  
  capture.output(Idents(function_object), file = "Cluster_IDs.txt")
  
  #After doing this clustering we perform the non-linear reduction of the data
  #using the same dimensions we already chose
  
  function_object <- RunUMAP(function_object, dims = selected_PC)
  
  #And we generate the DimPlot
  
  DPlot <- DimPlot(function_object, reduction = reduction_method)
  
  pdf(file = paste("DimPlot", reduction_method, "Neighbours.pdf", sep = "_"))
  print(DPlot)
  dev.off()
  
  ##Now we will use scCATCH to annotate the clusters
  
  library(scCATCH)
  
  plan(multicore, workers = cores)
  
  clu_markers <- findmarkergenes(function_object, species = organism) #Here we get the markers for each cluster
  clu_ann <- scCATCH(clu_markers$clu_markers, species = organism, tissue = selected_tissue) #And we get the annotation
  
  #Now we will print all the graphics
  #Starting with the violin plot for each of the genes found on findmarkergenes
  
  capture.output(paste("Marker genes are:", clu_markers$clu_markers, sep = "\n"))
  
  #Violinplot blablabla
  
  violin_counter = 1
  lista_genes = clu_markers$clu_markers$gene
  longitud = length(lista_genes)
  
  while(count <= longitud){
    print(VlnPlot(working_object, features = lista_genes[count]))
    violin_count <- violin_count +1
  }
  
  #Now we mark the Dimplot with the clusters
  
  new.cluster.ids <- as.list(clu_ann$cell_type)
  names(new.cluster.ids) <- levels(function_object)
  function_object <- RenameIdents(function_object, new.cluster.ids)
  
  Dimplotmarked <- DimPlot(function_object, reduction = reduction_method, label = TRUE, pt.size = 0.5) + NoLegend()
  
  pdf("UMAP_clustered.pdf")
  print(DimPlot(function_object, reduction = reduction_method, label = TRUE, pt.size = 0.5))
  dev.off()
  
  #And now we save the Robject so user can modify and inspect it on Rstudio
  
  saveRDS(function_object, file = "./Seurat_analysis.rds")
  
  return(function_object)
}

print("Clustering function... loaded!")

#After doing all this tasks we will have several .Rds files with the different steps
#and the date they were performed. We encourage the users to open any of the files and
#explore all the possiblities there


