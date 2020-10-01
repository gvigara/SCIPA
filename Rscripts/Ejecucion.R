project_name <- "pbmc3k"
input <- "/home/gvigast/Documentos/TFM/Pruebas/Prueba_Solo_R/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19" #esto vendr?a del script PERL
output <- "/home/gvigast/Documentos/TFM/Pruebas/Prueba_Solo_R/Analisis_Seurat"

date <- Sys.Date()
name <- paste("SeuratObject", project_name, date, sep="_")

#We gather the variables from the script


minimum_cells <- 3
minimum_features <- 200
minimum_subset_features <- 200
maximum_subset_features <- 2500
maximum_subset_mt <- 5

working_object <- setup_qc(input, output, project_name, minimum_cells,
                           minimum_features, minimum_subset_features, maximum_subset_features,
                           maximum_subset_mt)

saveRDS(working_object, file = paste(name, "QualityControl.rds", sep = "_"))

normalizacion <- "LogNormalize"
factor_normalizacion <- 10000
metodo_selecion <- "vst"
numero_caract <- 2000
caracter_top <- 10

working_object <- Normalization(working_object, normalizacion, factor_normalizacion,
                                  metodo_selecion, numero_caract, caracter_top)

saveRDS(working_object, file = paste(name, "NormalizedData.rds", sep = "_"))

dimensiones_texto <- 5
caract_texto <- 5
dimensiones_VDL <- 5
dimensiones_DHM <- 10
celulas_DHM <- 200

working_object <- PCA_and_dimensionality(working_object, dimensiones_texto, caract_texto,
                                         dimensiones_VDL, dimensiones_DHM, celulas_DHM)


saveRDS(working_object, file = paste(name, "PCA.rds", sep = "_"))

JackStraw_Replicates <- 100
JackStrawScore_dimensions <- 20
JackStrawPLot_dimensions <- 15

working_object <- dimensionality(working_object, JackStraw_Replicates,
                                 JackStrawScore_dimensions, JackStrawPLot_dimensions)

saveRDS(working_object, file = paste(name, "Dimensionalized.rds", sep = "_"))


eleccion_PC <- 10
resolucion_cluster <- 0.5
numero_celulas_ID <- 10
reduction_method <- "umap"

working_object <- clustering_and_NLR(working_object, eleccion_PC, resolucion_cluster,
                             numero_celulas_ID, reduction_method)

saveRDS(working_object, file = paste(name, "clustered.rds", sep = "_"))


#Now we generate the object, using the structure "SeuratObject_date"

