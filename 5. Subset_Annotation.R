# Subset datasets.integrated
# Recluster and reannotate (annotations.v1)
# Needs from script4: datasets.integrated with annotations.new.ordered

# Marker-gene modules
modules_merged <- list(
  Epiblast                  = c("OTX2", "DNMT3B", "NANOG", "SOX2", "POU5F1", "OCT4"),
  Primitive_Streak          = c("CDX2", "CDX1", "FGF8"),
  Anterior_Primitive_Streak = c("NODAL", "GSC", "MIXL1", "WNT3", "EOMES", "LEFTY1"),
  Node_Notochord            = c("SHH", "CHRD", "FOXA2", "GSC", "NOTO", "FOXJ1"),
  Caudal_Epiblast           = c("GBX2", "CDX1", "HOXA1", "FGF8", "CDX2"),
  Definitive_Endoderm       = c("SOX17", "FOXA2", "MIXL1", "HHEX", "GATA6", "CER1", "CHRD"),
  Nascent_Mesoderm          = c("MESP1", "TBXT", "CDX2", "CDX4", "FGF8", "WNT5A", "WNT8A", "WNT3A"),
  Cardiac_Mesoderm          = c("KDR", "HAND1", "MESP1", "BMP4"),
  Head_Mesoderm             = c("ALX1", "LIX1", "CYP26A1", "TCF21"),
  Presomitic_Mesoderm       = c("DLL1", "TBX6", "ALDH1A2", "MSGN1", "HES7"),
  Somitic_Mesoderm          = c("RIPPLY1", "RIPPLY2", "MESP2", "DLL3", "ALDH1A2", "FOXC2"),
  Dermomyotome              = c("DMRT2", "MEOX1", "PAX3", "MEOX2"),
  Sclerotome                = c("PAX1", "PAX9", "SOX9"),
  Connective_Mesenchyme     = c("PRRX2", "PDGFRB", "EYA1"),
  Myoblasts                 = c("MSC", "MYF5", "PAX7"),
  Myocytes_Myotubes         = c("MYOG", "MYMK", "MYOD1"),
  Myofibers                 = c("TECRL", "MYOT", "MYBPC3", "MYL2", "MYL1", "TNNT2", "ACTN2",
                                "MYH3", "MYH7", "MYH1", "MYL3"),
  Neural_Tube               = c("PAX6", "HES5", "NKX6-1", "NRCAM"),
  Neurons                   = c("SCN3A", "SYT4", "NEUROD1", "ISL2", "STMN4", "MAP2", "DCX", "STMN2")
)

# Subset object 
unique(datasets.integrated$annotations.new.ordered)
Idents(datasets.integrated) <- "annotations.new.ordered"
keep_levels <- c(
  "Neurons", "Epiblast/Neural plate", "Neural tube/Radial glia",
  "Caudal epiblast/Early neural tube", "Anterior primitive streak/Node",
  "Nascent mesoderm/Presomitic mesoderm", "Somitic mesoderm",
  "Dermomyotome", "Muscle")
atlas.subset <- subset(datasets.integrated, idents = keep_levels)
atlas.subset$old.annotations <- droplevels(factor(atlas.subset$annotations.new.ordered))
Idents(atlas.subset) <- "old.annotations"
unique(atlas.subset$old.annotations)

# Recluster object 
atlas.subset[["RNA"]] <- JoinLayers(atlas.subset[["RNA"]])
atlas.subset <- subset(atlas.subset, subset = nFeature_RNA > 2000)
DefaultAssay(atlas.subset) <- "integrated"
atlas.subset <- RunPCA(atlas.subset)
atlas.subset <- FindNeighbors(atlas.subset, dims = 1:30) 
atlas.subset <- FindClusters(atlas.subset, resolution = 1.8) 
atlas.subset$clusters_test <- Idents(atlas.subset)
atlas.subset <- RunUMAP(
  atlas.subset,
  dims = 1:32,
  n.neighbors = 60, 
  min.dist = 0.45, 
  spread = 0.8,
  return.model = TRUE,
  reduction.name = "umap_test"
)

# Module scores + bubble plot per cluster
mods <- lapply(modules_merged, intersect, rownames(atlas.subset[["RNA"]]))
mods <- mods[lengths(mods) >= 2]                       # drop modules with < 2 genes detected

for (nm in names(mods)) {
  atlas.subset <- AddModuleScore(atlas.subset, features = list(mods[[nm]]), name = "tmp",
                                 assay = "RNA", nbin = 24, ctrl = 100, seed = 1)
  atlas.subset[[paste0("MS_", nm)]] <- atlas.subset$tmp1
  atlas.subset$tmp1 <- NULL
}
ms_cols <- grep("^MS_", colnames(atlas.subset@meta.data), value = TRUE)

# Bubble plot of mean module score per group (by = a metadata column).
bubble <- function(by, colour_by_group = FALSE) {
  df <- atlas.subset@meta.data %>%
    mutate(grp = as.character(.data[[by]])) %>%
    group_by(grp) %>%
    summarise(across(all_of(ms_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_longer(-grp, names_to = "module", values_to = "mean_score") %>%
    mutate(module = sub("^MS_", "", module)) %>%
    group_by(grp) %>%
    mutate(z = as.numeric(scale(mean_score))) %>%
    ungroup()
  
  # groups in numeric order if they are numbers, else alphabetical
  lv <- unique(df$grp)
  lv <- if (!anyNA(suppressWarnings(as.numeric(lv)))) lv[order(as.numeric(lv))] else sort(lv)
  df$grp <- factor(df$grp, levels = lv)
  
  # modules ordered by the group in which they score highest
  tops <- df %>% group_by(grp) %>% slice_max(mean_score, n = 1, with_ties = FALSE) %>%
    arrange(grp) %>% pull(module)
  df$module <- factor(df$module, levels = c(unique(tops), setdiff(unique(df$module), tops)))
  
  p <- ggplot(df, aes(x = module, y = grp, size = mean_score)) +
    scale_size_area(max_size = 8) +
    coord_flip() +
    theme_classic() +
    theme(axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
  
  if (colour_by_group) {
    p + geom_point(aes(colour = grp)) + labs(x = "Module", y = by, size = "Mean score", colour = by)
  } else {
    p + geom_point(aes(colour = z)) + labs(x = "Module", y = by, size = "Mean score", colour = "Z-score")
  }
}

DimPlot(atlas.subset, reduction = "umap_test", group.by = "seurat_clusters", pt.size = 1) +
  bubble("seurat_clusters")


# Name the clusters: annotations.v1 
v1_by_cluster <- c(
  "0"  = "Epiblast 1",
  "1"  = "Nascent mesoderm/Cardiac mesoderm",
  "2"  = "Epiblast 1",
  "3"  = "Epiblast 2",
  "4"  = "Presomitic mesoderm",
  "5"  = "Primitive streak",
  "6"  = "Epiblast 2",
  "7"  = "Dermomyotome",
  "8"  = "Myoblasts",
  "9"  = "Nascent mesoderm",
  "10" = "Head mesoderm",
  "11" = "Neural tube/Radial glia",
  "12" = "Neural tube/Radial glia",
  "13" = "Nascent mesoderm/Cardiac mesoderm",
  "14" = "Presomitic mesoderm",
  "15" = "Primitive streak",
  "16" = "Epiblast 1",
  "17" = "Somitic mesoderm",
  "18" = "Node/Notochord",
  "19" = "Caudal epiblast",
  "20" = "Somitic connective mesenchyme",
  "21" = "Neurons",
  "22" = "Sclerotome",
  "23" = "Myofibers",
  "24" = "Neurons",
  "25" = "Somitic mesoderm",
  "26" = "Anterior primitive streak",
  "27" = "Definitive endoderm",
  "28" = "Myocytes/Myotubes",
  "29" = "Myofibers"
)

v1_levels <- c("Epiblast 1", "Nascent mesoderm/Cardiac mesoderm", "Head mesoderm",
               "Neural tube/Radial glia", "Presomitic mesoderm", "Primitive streak",
               "Somitic mesoderm", "Node/Notochord", "Caudal epiblast",
               "Somitic connective mesenchyme", "Neurons", "Sclerotome", "Myofibers",
               "Anterior primitive streak", "Definitive endoderm", "Myocytes/Myotubes",
               "Epiblast 2", "Dermomyotome", "Myoblasts", "Nascent mesoderm")

atlas.subset$annotations.v1 <- factor(unname(v1_by_cluster[as.character(atlas.subset$seurat_clusters)]),
                                      levels = v1_levels)
Idents(atlas.subset) <- "annotations.v1"
table(atlas.subset$annotations.v1)
DimPlot(atlas.subset, reduction = "umap_test", group.by = "annotations.v1", label = TRUE, repel = TRUE)

