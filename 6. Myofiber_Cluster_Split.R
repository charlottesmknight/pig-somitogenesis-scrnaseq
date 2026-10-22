# Separate cardiac muscle from skeletal myofibers in cluster 23 (annotations.v2)
# Needs from script 5 atlas.subset with seurat_clusters and annotations.v1

DefaultAssay(atlas.subset) <- "RNA"
is23 <- atlas.subset$seurat_clusters == "23"
is29 <- atlas.subset$seurat_clusters == "29"

# Candidate gene sets, filtered against cluster 29 (cardiac genes absent, skeletal genes detected)
cardiac_candidates  <- c("NKX2-5", "GATA4", "TBX5", "TBX20", "HAND2",
                         "MYBPC3", "RYR2", "TNNI3", "PLN", "MYL7", "NPPA", "MYH6")
skeletal_candidates <- c("MYH3", "MYL1", "MYLPF", "CKM", "TNNT3", "TNNI2",
                         "TNNC2", "ATP2A1", "MYBPC2", "MYOG", "MYF6")

expr <- GetAssayData(atlas.subset, assay = "RNA", layer = "data")
cardiac_present  <- intersect(cardiac_candidates, rownames(expr))
skeletal_present <- intersect(skeletal_candidates, rownames(expr))
pct29 <- 100 * Matrix::rowMeans(expr[c(cardiac_present, skeletal_present), is29] > 0)
cardiac_use  <- cardiac_present[pct29[cardiac_present] <= 5]
skeletal_use <- skeletal_present[pct29[skeletal_present] >= 10]

# Module scores and split of cluster 23
atlas.subset <- AddModuleScore(atlas.subset, features = list(cardiac_use),  name = "Cardiac_core",  seed = 1)
atlas.subset <- AddModuleScore(atlas.subset, features = list(skeletal_use), name = "Skeletal_core", seed = 1)

cardiac_cut  <- 0.3
skeletal_cut <- 0.6
is_cardiac <- is23 & atlas.subset$Cardiac_core1 > cardiac_cut & atlas.subset$Skeletal_core1 < skeletal_cut

# annotations.v2: cardiac cells become separate cluster
v2 <- as.character(atlas.subset$annotations.v1)
v2[is_cardiac] <- "Cardiac muscle"
v1_levels <- levels(atlas.subset$annotations.v1)
atlas.subset$annotations.v2 <- factor(v2, levels = append(v1_levels, "Cardiac muscle", after = which(v1_levels == "Myofibers")))
table(atlas.subset$annotations.v2)

Idents(atlas.subset) <- "annotations.v2"
saveRDS(atlas.subset, "atlas.subset.final.v2.rds")

