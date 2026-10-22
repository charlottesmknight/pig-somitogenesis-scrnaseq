# === Integration === 
# Needs Pig_atlas_Simpson_2024.RDS (Pig gastrulation atlas, Simpson et al., 2024)
# Needs Myo_data_Cai_2023.RDS (Pig myogenesis data, Cai et al., 2023)
# Produces: datasets.integrated (integrated, clustered object), pig_gastrulation_dataset and pig_myogenesis_data (rebuilt, filtered source objects)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)     
  library(tibble)
  library(ggplot2)
  library(pheatmap)  
})

pig_myogenesis_data <- readRDS("Myo_data_Cai_2023.RDS")
pig_gastrulation_dataset <- readRDS("Pig_atlas_Simpson_2024.RDS")

pig_myogenesis_data$dataset <- "pig_myogenesis_data"
pig_gastrulation_dataset$dataset <- "pig_gastrulation_dataset"

pig_myo_data <- GetAssayData(pig_myogenesis_data, assay = "RNA")
pig_gast_data <- GetAssayData(pig_gastrulation_dataset, assay = "RNA")

# === Restrict atlas to genes shared with myogenesis dataset === 
common_genes <- intersect(rownames(pig_gast_data), rownames(pig_myo_data))
pig_gast_data <- pig_gast_data[common_genes, ]

# === Run standard workflow === 
DefaultAssay(pig_gastrulation_dataset) <- "RNA" 
pig_gast_meta_data <- pig_gastrulation_dataset@meta.data
pig_gastrulation_dataset <- CreateSeuratObject(counts = pig_gast_data, project = "pig_gastrulation_dataset",
                                               meta.data = pig_gast_meta_data, min.cells = 0, min.features = 1750)
pig_gastrulation_dataset <- NormalizeData(pig_gastrulation_dataset, verbose = FALSE)
pig_gastrulation_dataset <- FindVariableFeatures(pig_gastrulation_dataset, verbose = FALSE)
pig_gastrulation_dataset <- ScaleData(pig_gastrulation_dataset, features = rownames(pig_gastrulation_dataset), verbose = FALSE)
pig_gastrulation_dataset <- RunPCA(pig_gastrulation_dataset, verbose = FALSE, npcs = 100)

DefaultAssay(pig_myogenesis_data) <- "RNA" 
pig_myo_meta_data <- pig_myogenesis_data@meta.data
pig_myogenesis_data <- CreateSeuratObject(counts = pig_myo_data, project = "pig_myogenesis_data", 
                                          meta.data = pig_myo_meta_data, min.cells = 0, min.features = 1250)
pig_myogenesis_data <- NormalizeData(pig_myogenesis_data, verbose = FALSE)
pig_myogenesis_data <- FindVariableFeatures(pig_myogenesis_data, verbose = FALSE)
pig_myogenesis_data <- ScaleData(pig_myogenesis_data, features = rownames(pig_myogenesis_data), verbose = FALSE)
pig_myogenesis_data <- RunPCA(pig_myogenesis_data, verbose = FALSE, npcs = 100)


# === Label transfer, using gastrulation dataset as reference ===
transfer.anchors <- FindTransferAnchors(reference = pig_gastrulation_dataset, query = pig_myogenesis_data, 
                                        dims = 1:75, reference.reduction = "pca", n.trees = 600)
cellsubtype.predictions <- TransferData(anchorset = transfer.anchors, refdata = pig_gastrulation_dataset$Cellsubclusters, 
                                        dims = 1:75, n.trees = 600)
pig_myogenesis_data <- AddMetaData(pig_myogenesis_data, metadata = cellsubtype.predictions)
pig_myogenesis_data$predicted_Pig_atlas_cellsubtype <- pig_myogenesis_data$predicted.id
pig_myogenesis_data <- FindNeighbors(pig_myogenesis_data, reduction = "pca", dims = 1:25, return.neighbor = T)
pig_myogenesis_data <- RunUMAP(pig_myogenesis_data, reduction = "pca", dims = 1:25, n.neighbors = 30L) #30 pcs
DimPlot(pig_myogenesis_data, group.by = "predicted_Pig_atlas_cellsubtype")


# === Run Integration ===
object_list <- c(SplitObject(pig_myogenesis_data, split.by = "orig.ident"), SplitObject(pig_gastrulation_dataset, split.by = "Sample"))
object_list <- lapply(X = object_list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
  x <- ScaleData(x, features = rownames(x), verbose = FALSE)
  x <- RunPCA(x, verbose = FALSE, npcs = 100)
})

features <- SelectIntegrationFeatures(object.list = object_list, nfeatures = 2000) 
anchors <- FindIntegrationAnchors(object.list = object_list, anchor.features = features, reference = NULL, reduction = "rpca", dims = 1:75,  n.trees = 600)
datasets.integrated <- IntegrateData(anchorset = anchors)

# === Run standard workflow on integrated data ===
datasets.integrated[["RNA"]] <- JoinLayers(datasets.integrated[["RNA"]])
DefaultAssay(datasets.integrated) <- "integrated"
datasets.integrated <- subset(datasets.integrated, subset = nFeature_RNA > 2000)
datasets.integrated <- ScaleData(datasets.integrated)
datasets.integrated <- FindVariableFeatures(datasets.integrated)
datasets.integrated <- RunPCA(datasets.integrated, npcs = 50, verbose = FALSE)
datasets.integrated <- FindNeighbors(datasets.integrated, reduction = "pca", dims = 1:30)
datasets.integrated <- FindClusters(datasets.integrated, resolution = 0.5)
datasets.integrated <- RunUMAP(datasets.integrated, reduction = "pca", dims = 1:30)
datasets.integrated$clusters_after_integration <- datasets.integrated$seurat_clusters

# === Add the predicted cell types in the integrated dataset === 
datasets.integrated$Gastrulation_celltype_predictions <- ifelse(is.na(datasets.integrated$Cellsubclusters), 
                                                                datasets.integrated$predicted_Pig_atlas_cellsubtype,
                                                                datasets.integrated$Cellsubclusters)
DimPlot(datasets.integrated, group.by = "Gastrulation_celltype_predictions")



