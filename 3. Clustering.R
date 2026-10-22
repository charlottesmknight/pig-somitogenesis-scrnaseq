# Cluster with res 0.8, prune.SNN = 0.1
# Final UMAP parameters reduction.new.params
# Needs from script 2: datasets.integrated

datasets.integrated[["RNA"]] <- JoinLayers(datasets.integrated[["RNA"]])
DefaultAssay(datasets.integrated) <- "integrated"
datasets.integrated <- subset(datasets.integrated, subset = nFeature_RNA > 2000)
datasets.integrated <- ScaleData(datasets.integrated)
datasets.integrated <- FindVariableFeatures(datasets.integrated)
datasets.integrated <- RunPCA(datasets.integrated, npcs = 75, verbose = FALSE)
datasets.integrated <- FindNeighbors(datasets.integrated, reduction="pca", dims=1:25, prune.SNN=1/10)
datasets.integrated <- FindClusters(datasets.integrated, resolution=0.8)
datasets.integrated$clusters.new.params <- Idents(datasets.integrated)
datasets.integrated <- RunUMAP(datasets.integrated,dims = 1:30, 
                               min.dist=0.6, seed.use=123, 
                               reduction.name="reduction.new.params")

DimPlot(datasets.integrated, group.by = "clusters.new.params", reduction = "reduction.new.params")
length(unique(datasets.integrated$clusters.new.params))
