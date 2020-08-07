#Seurat R script

#Get the working directory
#Here we should pass the workign directory variable from PERL
#The software should ask for where they want to store the data

setwd("/home/gvigast/Documentos/Pruebas_Scritps_Seurat/pbmc_1k_v3_fastqs/Seurat_analysis")

#Load the PBMC dataset with the data from Cellranger
#We should add the directory to the script from perl
pbmc.data <- Read10X(data.dir = "/home/gvigast/Documentos/Pruebas_Scritps_Seurat/pbmc_1k_v3_fastqs/pbmc1k/outs/filtered_feature_bc_matrix");

#Initialize the Raw Seurat Object

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc1k", min.cells = 3, min.features = 200);

#After the object has been created we initiate the standard pre-processing
#in which we obtain the mitochondrial DNA

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-");

#Now we export the graphs with different visualization

#We create a violin plot of the data and export it to pdf

pdf("violin_plot_mitocount_pbmc1k.pdf"); #This generates the file
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3);
dev.off(); #This closes the graph

#Also, we create a scatter plot of the features

scatter_plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt");
scatter_plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 ="nFeature_RNA");

#After we generate the plots we represent them separately and together

pdf("scatter_plot_features_pbmc1k.pdf");
scatter_plot1;
scatter_plot2;
scatter_plot1 + scatter_plot2;
dev.off();

#And now we determine the working subset with the features selected
#by the user. Here we will use the ones by defect, which are the following
#But we have to substitute the values for those determined by the user

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt <5);

# THIS SCRIPT ENDS HERE 
