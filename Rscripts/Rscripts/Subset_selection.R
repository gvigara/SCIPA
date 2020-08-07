## Data Normalization and feature selection

#First we want to normalize the data, by using
#the function NormalizeData. The usual method uses LogNormalize,
#althoug there are several methods avaiable

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

#We identify the top10 highly variable genes

top10 <- head(VariableFeatures(pbmc), 10);

#After this, we obtain the variable features. By defect, the selection method
#is vst although there are others (notes on the program) and the default
#number of features is 2000

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

#Now we plot the variable features in a picture (PDF file)

plot1VF <- VariableFeaturePlot(pbmc);
plot2VF <- LabelPoints(plot = plot1VF, points = top10, repel = TRUE);
pdf(file="Variablefeatures_plot.pdf");
plot1VF;
plot2VF;
dev.off();

#END OF FILE
