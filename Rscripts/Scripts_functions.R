## SCRIPT: ANALYSIS OF scRNAseq Data

libraries <- c("Seurat", "dplyr", "patchwork", "future")
sapply(libraries, require, character.only=TRUE)

setup_qc <- function(input_directory, output_directory, project_name, min_cells, min_features, 
                     min_nFeature_subset = "." , max_nFeature_subset = ".", max_percent_mt =".", organism){
  
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
  
  #If it is human the pattern for mitochondrial DNA is ^MT-, it mouse is ^mt-
  
  if (organism == "Human"){
  project_name[["percent.mt"]] <- PercentageFeatureSet(project_name, pattern = "^MT-")
  }
  else if (organism == "human"){
  project_name[["percent.mt"]] <- PercentageFeatureSet(project_name, pattern = "^MT-")
  }
  else if (organism == "Mouse"){
    project_name[["percent.mt"]] <- PercentageFeatureSet(project_name, pattern = "^mt-")	
  }
  else if (organism == "mouse"){
    project_name[["percent.mt"]] <- PercentageFeatureSet(project_name, pattern = "^mt-")	
  }
  else{
  print("Organism selected is not recognizable, please check config.ini")
  exit()
  }
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
  
  if(min_nFeature_subset == "." && max_nFeature_subset == "." && max_percent_mt == "."){
  minimum_subset_features = min(project_name[["nFeature_RNA"]]$nFeature_RNA)
  #As maximum subset we will use the mean plus the intercuartilic range of the features
  maximum_subset_features = (1.5 * IQR(project_name[["nFeature_RNA"]]$nFeature_RNA) + quantile(project_name[["nFeature_RNA"]]$nFeature_RNA, 0.75)
  #As mitochondrial value we will use the values under the percentile 99, gathering most of the values there
  maximum_mitochondral_dna_percentage = (1.5 * IQR(project_name[["percent.mt"]]$percent.mt) + quantile(project_name[["percent.mt"]]$percent.mt, 0.75)
  print("The user didn't provide any parameters for subset generation, using the calculated ones")
  print(paste("Minimum subset features:", minimum_subset_features, sep = " "))
  print(paste("Maximum_subset_features:", maximum_subset_features, sep = " "))
  print(paste("Maximum_mitochondrial_dna_percentage:", maximum_mitochondral_dna_percentage, sep = " "))
  project_name <- subset(project_name, subset = nFeature_RNA > minimum_subset_features & nFeature_RNA < maximum_subset_features & percent.mt < maximum_mitochondral_dna_percentage)
  } 
  else if(min_nFeature_subset != "." && max_nFeature_subset == "." && max_percent_mt == "."){
  minimum_subset_features = min_nFeature_subset
  #As maximum subset we will use the mean plus the intercuartilic range of the features
  maximum_subset_features = IQR(project_name[["nFeature_RNA"]]$nFeature_RNA) + mean(project_name[["nFeature_RNA"]]$nFeature_RNA)
  #As mitochondrial value we will use the values under the percentile 99, gathering most of the values there
  maximum_mitochondral_dna_percentage = IQR(project_name[["percent.mt"]]$percent.mt) + mean(project_name[["percent.mt"]]$percent.mt)
  print("The user didn't provide any parameters for subset generation, using the calculated ones")
  print(paste("Minimum subset features:", minimum_subset_features, sep = " "))
  print(paste("Maximum_subset_features:", maximum_subset_features, sep = " "))
  print(paste("Maximum_mitochondrial_dna_percentage:", maximum_mitochondral_dna_percentage, sep = " "))
  project_name <- subset(project_name, subset = nFeature_RNA > minimum_subset_features & nFeature_RNA < maximum_subset_features & percent.mt < maximum_mitochondral_dna_percentage)
  } 
  else if(min_nFeature_subset == "." && max_nFeature_subset != "." && max_percent_mt == "."){
  minimum_subset_features = min(project_name[["nFeature_RNA"]]$nFeature_RNA)
  #As maximum subset we will use the mean plus the intercuartilic range of the features
  maximum_subset_features = max_nFeature_subset
  #As mitochondrial value we will use the values under the percentile 99, gathering most of the values there
  maximum_mitochondral_dna_percentage = IQR(project_name[["percent.mt"]]$percent.mt) + mean(project_name[["percent.mt"]]$percent.mt)
  print("The user didn't provide any parameters for subset generation, using the calculated ones")
  print(paste("Minimum subset features:", minimum_subset_features, sep = " "))
  print(paste("Maximum_subset_features:", maximum_subset_features, sep = " "))
  print(paste("Maximum_mitochondrial_dna_percentage:", maximum_mitochondral_dna_percentage, sep = " "))
  project_name <- subset(project_name, subset = nFeature_RNA > minimum_subset_features & nFeature_RNA < maximum_subset_features & percent.mt < maximum_mitochondral_dna_percentage)
  }
  else if(min_nFeature_subset == "." && max_nFeature_subset == "." && max_percent_mt != "."){
  minimum_subset_features = min(project_name[["nFeature_RNA"]]$nFeature_RNA)
  #As maximum subset we will use the mean plus the intercuartilic range of the features
  maximum_subset_features = IQR(project_name[["nFeature_RNA"]]$nFeature_RNA) + mean(project_name[["nFeature_RNA"]]$nFeature_RNA)
  #As mitochondrial value we will use the values under the percentile 99, gathering most of the values there
  maximum_mitochondral_dna_percentage = max_percent_mt
  print("The user didn't provide any parameters for subset generation, using the calculated ones")
  print(paste("Minimum subset features:", minimum_subset_features, sep = " "))
  print(paste("Maximum_subset_features:", maximum_subset_features, sep = " "))
  print(paste("Maximum_mitochondrial_dna_percentage:", maximum_mitochondral_dna_percentage, sep = " "))
  project_name <- subset(project_name, subset = nFeature_RNA > minimum_subset_features & nFeature_RNA < maximum_subset_features & percent.mt < maximum_mitochondral_dna_percentage)
  } 
  else if(min_nFeature_subset != "." && max_nFeature_subset != "." && max_percent_mt == "."){
  minimum_subset_features = min_nFeature_subset
  #As maximum subset we will use the mean plus the intercuartilic range of the features
  maximum_subset_features = max_nFeature_subset
  #As mitochondrial value we will use the values under the percentile 99, gathering most of the values there
  maximum_mitochondral_dna_percentage = IQR(project_name[["percent.mt"]]$percent.mt) + mean(project_name[["percent.mt"]]$percent.mt)
  print("The user didn't provide any parameters for subset generation, using the calculated ones")
  print(paste("Minimum subset features:", minimum_subset_features, sep = " "))
  print(paste("Maximum_subset_features:", maximum_subset_features, sep = " "))
  print(paste("Maximum_mitochondrial_dna_percentage:", maximum_mitochondral_dna_percentage, sep = " "))
  project_name <- subset(project_name, subset = nFeature_RNA > minimum_subset_features & nFeature_RNA < maximum_subset_features & percent.mt < maximum_mitochondral_dna_percentage)
  } 
  else if(min_nFeature_subset != "." && max_nFeature_subset == "." && max_percent_mt != "."){
  minimum_subset_features = min_nFeature_subset
  #As maximum subset we will use the mean plus the intercuartilic range of the features
  maximum_subset_features = IQR(project_name[["nFeature_RNA"]]$nFeature_RNA) + mean(project_name[["nFeature_RNA"]]$nFeature_RNA)
  #As mitochondrial value we will use the values under the percentile 99, gathering most of the values there
  maximum_mitochondral_dna_percentage = max_percent_mt
  print("The user didn't provide any parameters for subset generation, using the calculated ones")
  print(paste("Minimum subset features:", minimum_subset_features, sep = " "))
  print(paste("Maximum_subset_features:", maximum_subset_features, sep = " "))
  print(paste("Maximum_mitochondrial_dna_percentage:", maximum_mitochondral_dna_percentage, sep = " "))
  project_name <- subset(project_name, subset = nFeature_RNA > minimum_subset_features & nFeature_RNA < maximum_subset_features & percent.mt < maximum_mitochondral_dna_percentage)
  }
  else if(min_nFeature_subset == "." && max_nFeature_subset != "." && max_percent_mt != "."){
  minimum_subset_features = min(project_name[["nFeature_RNA"]]$nFeature_RNA)
  #As maximum subset we will use the mean plus the intercuartilic range of the features
  maximum_subset_features = max_nFeature_subset
  #As mitochondrial value we will use the values under the percentile 99, gathering most of the values there
  maximum_mitochondral_dna_percentage = max_percent_mt
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
  #saveRDS(project_name, file = 'Analysis_seurat_object.rds')
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
  jpeg(filename = "Seurat_VarFeaturePlot.jpg", units ="px", width = 2048, height = 1024)  #Opens the picture
  print(NFeatures_plot1 + NFeatures_plot2)
  dev.off()
  
  
  return(function_object)
  
}

print("Normalization function... loaded!")

## Third function will be in charge of PCA and determination of the 
## dimensionality of the dataset

PCA <- function(parse_object, dimensions_on_text,
                                   features_on_text, VDL_dims, DHM_dims,
                                   DHM_cells, cores){
  
  function_object <- parse_object
  plan(multicore, workers = cores)
  
  print("Starting PCA adjustment:")
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
                       nfeatures = features_on_text), file = "Genes_on_PCAs.txt")
  
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
  
  pdf(file = "UMAP_nocluster.pdf")
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
                               reduction_method, cores){
  
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
  
  #capture.output(
  #paste(print("Selected_PCs_are:"),
  #print(selected_PC), sep = ""), file = "Selected_dimensions.txt")

  
  #We first find the neighbors of the clusters. The dimensions are the ones 
  #described by the user previously based on the information obtained by JackStraw
  
  function_object <- FindNeighbors(function_object, dims = selected_PC)
  
  #And after this we obtain the clusters. The value needed is the resolution,
  #that will be provided with the variable cluster_resolution
  
  function_object <- FindClusters(function_object, resolution = cluster_resolution)
  
  #With this information we can see the IDs of the first N cells asked by the user
  #in the cluster
  
  options(max.print= 99999999) #With this we avoid the limit for printing
  
  capture.output(Idents(function_object), file = "Cluster_cellIDs.txt")
  
  #After doing this clustering we perform the non-linear reduction of the data
  #using the same dimensions we already chose
  
  function_object <- RunUMAP(function_object, dims = selected_PC)
  
  #And we generate the DimPlot
  
  DPlot <- DimPlot(function_object, reduction = reduction_method)
  
  pdf(file = "UMAP_clustered_unidentified.pdf"))
  print(DPlot)
  dev.off()
  
  return(function_object)
}

print("Clustering function... loaded!")

scCATCH_annotation <- function(parse_object, organism, selected_tissue, plotted_genes_cluster, heatmap_genes, cores){

  ##Now we will use scCATCH to annotate the clusters
  
  cat("\nInitiating scCATCH auto-typing of clusters: \n")
  print(paste("Organism for annotation:", organism, sep = " "))
  print(paste("Selected tissue for annotation:", selected_tissue, sep = " "))
  
  library(scCATCH)
  
  plan(multicore, workers = cores)
  
  function_object <- parse_object
  
  clu_markers <- findmarkergenes(function_object, species = organism) #Here we get the markers for each cluster
  
  saveRDS(clu_markers, file = "scCATCH_cluster_markers.rds")
  
  #we capture the marker genes for all the clusters:
  write.csv(clu_markers$clu_markers, 'scCATCH_cluster_markers.csv')
  #capture.output(paste("Marker genes are:", clu_markers$clu_markers, sep = "\n"), file = "cluster_markers.txt")
  
  #We print the violin plot 
  
  violin_counter = 1
  lista_genes_violinplot = clu_markers$clu_markers %>% group_by(cluster) %>% top_n(n = plotted_genes_cluster, wt = avg_logfc) #we select the top n genes for each cluster
  write.csv(lista_genes_violinplot, 'scCATCH_ViolinPlot.csv') #write the list of cluster markers to a file 
  lista_genes = lista_genes_violinplot$gene
  
  longitud = length(as.list(lista_genes))
  
  #Print the violin plot for each marker gene
  pdf(file = "scCATCH_ViolinPlot.pdf")
  while(violin_counter <= longitud){
    print(VlnPlot(function_object, features = lista_genes[violin_counter]))
    violin_counter <- violin_counter +1
  }
  dev.off()
  
  #Print the feature plot for each marker gene
  violin_counter = 1 #we reset the counter
  pdf(file = "scCATCH_FeaturePlot.pdf")
  while(violin_counter <= longitud){
    print(FeaturePlot(function_object, features = lista_genes[violin_counter]))
    violin_counter <- violin_counter +1
  }
  dev.off()
  
  #we print the heatmap genes for the clustering
  
  topN = clu_markers$clu_markers %>% group_by(cluster) %>% top_n(n = heatmap_genes, wt = avg_logfc) #we select the top n genes for each cluster
  write.csv(topN, 'scCATCH_Heatmap.csv') #write the list of cluster markers to a file
  lista_genes_heatmap = lista_genes_violinplot$gene
  pdf(file = "scCATCH_Heatmap.pdf")
  print(DoHeatmap(function_object, features = lista_genes_heatmap) + NoLegend())
  dev.off()
  
  #And finally we annotate the clusters
  clu_ann <- scCATCH(clu_markers$clu_markers, species = organism, tissue = selected_tissue) #And we get the annotation
  
  return(clu_ann)
}

print("scCATCH typing function... loaded!")


cluster_marking <- function(parse_object, scCATCH_data, reduction_method){
  #Now we implement the function where unmatched cell types are typed as Unknown
  
  function_object <- parse_object
  clu_ann <- scCATCH_data #this is the annotation by scCATCH done in the previous function
  
  tmp1 <- data.frame(cluster = levels(Idents(function_object)))
  tmp <- merge(tmp1, clu_ann, by = 'cluster', all = T)
  tmp$cell_type[which(is.na(tmp$cell_type))] <- "Unclassified"

  #And we assign the clusters
  
  new.cluster.ids <- tmp$cell_type
  names(new.cluster.ids) <- levels(function_object)
  function_object <- RenameIdents(function_object, new.cluster.ids)
  
  #Now we mark the Dimplot with the clusters
  
  Dplot2 <- DimPlot(function_object, reduction = reduction_method, label = TRUE, pt.size = 0.5)
  
  pdf(file = "UMAP_clustered.pdf")  #Opens the picture
  print(Dplot2)
  dev.off()
  
  Dplot3 <- DimPlot(function_object, reduction = reduction_method, label = TRUE, pt.size = 0.5) + NoLegend()
  
  pdf(file = "UMAP_clustered_nolegend.pdf")
  print(Dplot3)
  dev.off()
  
  Dplot4 <- DimPlot(function_object, reduction = reduction_method, label = FALSE, pt.size = 0.5)
  
  pdf(file = "UMAP_clustered_nolabel.pdf")
  print(Dplot4)
  dev.off()
  
  saveRDS(function_object, file = "Seurat_analysis_final.rds")
  
  cat("RDS file saved. Function terminated\n")
  
  return(function_object)
  
}

print("scCATCH annotation function... loaded!")


#After doing all this tasks we will have several .Rds files with the different steps
#and the date they were performed. We encourage the users to open any of the files and
#explore all the possiblities there


