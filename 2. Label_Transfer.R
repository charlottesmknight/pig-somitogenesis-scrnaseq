# === Recluster ===
# Needs, from script 1: datasets.integrated, pig_gastrulation_dataset, pig_myogenesis_data
# Needs Marioni_Atlas_pig_names_V5.RDS (mouse gastrulation atlas, Pijuan-Sala et al., 2019) with mouse gene names converted to pig gene names

DefaultAssay(datasets.integrated) <- "integrated"
datasets.integrated <- subset(datasets.integrated, subset = nFeature_RNA > 2000)
datasets.integrated <- FindVariableFeatures(datasets.integrated)
datasets.integrated <- ScaleData(datasets.integrated, features = rownames(datasets.integrated))
datasets.integrated <- RunPCA(datasets.integrated, npcs = 50, verbose = FALSE)
datasets.integrated <- RunUMAP(datasets.integrated, reduction = "pca", dims = 1:30, n.neighbors = 30L) #30 pcs
datasets.integrated <- FindNeighbors(datasets.integrated, dims = 1:30, reduction = "pca", k.param = 30)
datasets.integrated <- FindClusters(datasets.integrated, resolution = 1.2, group.singletons = FALSE)
datasets.integrated$seurat.clusters.post.integration <- datasets.integrated$integrated_snn_res.1.2
datasets.integrated <- SetIdent(datasets.integrated, value = "seurat.clusters.post.integration")


# === Label transfer, using gastrulation dataset as reference === 

DefaultAssay(pig_gastrulation_dataset) <- "RNA" 
pig_gastrulation_dataset <- NormalizeData(pig_gastrulation_dataset, verbose = FALSE)
pig_gastrulation_dataset <- FindVariableFeatures(pig_gastrulation_dataset, verbose = FALSE)
pig_gastrulation_dataset <- ScaleData(pig_gastrulation_dataset, features = rownames(pig_gastrulation_dataset), verbose = FALSE)
pig_gastrulation_dataset <- RunPCA(pig_gastrulation_dataset, verbose = FALSE, npcs = 100)

DefaultAssay(datasets.integrated) <- "RNA" 
datasets.integrated <- NormalizeData(datasets.integrated, verbose = FALSE)
datasets.integrated <- FindVariableFeatures(datasets.integrated, verbose = FALSE)
datasets.integrated <- ScaleData(datasets.integrated, features = rownames(datasets.integrated), verbose = FALSE)
datasets.integrated <- RunPCA(datasets.integrated, verbose = FALSE, npcs = 100)

transfer.anchors.gast.as.ref <- FindTransferAnchors(reference = pig_gastrulation_dataset, 
                                        query = datasets.integrated, 
                                        dims = 1:75, reference.reduction = "pca", n.trees = 600)
predictions.gast.as.ref <- TransferData(anchorset = transfer.anchors.gast.as.ref, 
                                        refdata = pig_gastrulation_dataset$Cellsubclusters, 
                                        dims = 1:75, n.trees = 600)
datasets.integrated <- AddMetaData(datasets.integrated, metadata = predictions.gast.as.ref)
datasets.integrated$predictions.gast <- datasets.integrated$predicted.id



# === Label transfer, using myogenesis dataset as reference === 

DefaultAssay(pig_myogenesis_data) <- "RNA" 
pig_myogenesis_data <- NormalizeData(pig_myogenesis_data, verbose = FALSE)
pig_myogenesis_data <- FindVariableFeatures(pig_myogenesis_data, verbose = FALSE)
pig_myogenesis_data <- ScaleData(pig_myogenesis_data, features = rownames(pig_myogenesis_data), verbose = FALSE)
pig_myogenesis_data <- RunPCA(pig_myogenesis_data, verbose = FALSE, npcs = 100)

transfer.anchors.myo.as.ref <- FindTransferAnchors(reference = pig_myogenesis_data, 
                                                    query = datasets.integrated, 
                                                    dims = 1:75, reference.reduction = "pca", n.trees = 600)
predictions.myo.as.ref <- TransferData(anchorset = transfer.anchors.myo.as.ref, 
                                        refdata = pig_myogenesis_data$Myogenesis_celltypes, 
                                        dims = 1:75, n.trees = 600)
datasets.integrated <- AddMetaData(datasets.integrated, metadata = predictions.myo.as.ref)
datasets.integrated$predictions.myo <- datasets.integrated$predicted.id


# === Label transfer, using Marioni dataset as reference === 

mouse_atlas_marioni <- readRDS("Marioni_Atlas_pig_names_V5.RDS")
DefaultAssay(mouse_atlas_marioni) <- "RNA" 
mouse_atlas_marioni <- NormalizeData(mouse_atlas_marioni, verbose = FALSE)
mouse_atlas_marioni <- FindVariableFeatures(mouse_atlas_marioni, verbose = FALSE)
mouse_atlas_marioni <- ScaleData(mouse_atlas_marioni, features = rownames(mouse_atlas_marioni), verbose = FALSE)
mouse_atlas_marioni <- RunPCA(mouse_atlas_marioni, verbose = FALSE, npcs = 100)

datasets.integrated <- SetIdent(datasets.integrated, value = "seurat.clusters.post.integration")
mouse_atlas_marioni <- SetIdent(mouse_atlas_marioni, value = "celltype_extended_atlas")

transfer.anchors.mouse.as.ref <- FindTransferAnchors(reference = mouse_atlas_marioni, 
                                                   query = datasets.integrated, 
                                                   dims = 1:50, reference.reduction = "pca", n.trees = 600)
predictions.mouse.as.ref <- TransferData(anchorset = transfer.anchors.mouse.as.ref, 
                                       refdata = mouse_atlas_marioni$celltype_extended_atlas, 
                                       dims = 1:50, n.trees = 600)
datasets.integrated <- AddMetaData(datasets.integrated, metadata = predictions.mouse.as.ref)
datasets.integrated$predictions.mouse <- datasets.integrated$predicted.id


# === Visualise label transfer === 

DefaultAssay(datasets.integrated) <- "integrated"
DimPlot(datasets.integrated, group.by = "predictions.myo") + theme(legend.position = "bottom") + coord_fixed()
DimPlot(datasets.integrated, group.by = "predictions.gast") + theme(legend.position = "bottom") + coord_fixed()
DimPlot(datasets.integrated, group.by = "predictions.mouse") + theme(legend.position = "bottom") + coord_fixed()
DimPlot(datasets.integrated, group.by = "seurat.clusters.post.integration") + theme(legend.position = "bottom") + coord_fixed()


# === assigning transferred labels === 

threshold <- 0.75
myo_ds    <- "pig_myogenesis_data"

prop_myo <- tapply(
  datasets.integrated@meta.data[["dataset"]] == myo_ds,                   
  datasets.integrated@meta.data[["seurat.clusters.post.integration"]],      
  mean, na.rm = TRUE                                                        
)
clusters_use_myo <- names(prop_myo)[prop_myo >= threshold]
use_myo <- datasets.integrated@meta.data[["seurat.clusters.post.integration"]] %in% clusters_use_myo
datasets.integrated[["predictions.merged"]] <- ifelse(
  use_myo,
  datasets.integrated@meta.data[["predictions.myo"]],
  datasets.integrated@meta.data[["predictions.gast"]]
)
DimPlot(datasets.integrated, group.by = "predictions.merged", label = TRUE, repel = TRUE) + NoLegend()


