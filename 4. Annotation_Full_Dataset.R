# Annotation on res 0.8 clusters (annotations.new.ordered)
# Needs from script 3: datasets.integrated with clusters.new.params, reduction.new.params and predictions.gast

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)    
  library(ggplot2)
  library(patchwork) 
  library(pheatmap)
})

datasets.integrated$annotations.new <- datasets.integrated$clusters.new.params
datasets.integrated <- SetIdent(datasets.integrated, value = "annotations.new")
unique(datasets.integrated$clusters.new.params)

# Cluster markers
markers_res0.8 <- FindAllMarkers(datasets.integrated, assay = "RNA", only.pos = TRUE,
                                 min.pct = 0.25, logfc.threshold = 0.25)

datasets.integrated <- RenameIdents(datasets.integrated, 
                                    "0"  = "Mesenchyme 3", 
                                    "1"  = "Nascent mesoderm/Presomitic mesoderm",
                                    "2"  = "Mesenchyme 1",
                                    "3"  = "Epiblast/Neural plate", 
                                    "4"  = "Neural tube/Radial glia", 
                                    "5"  = "Lateral plate mesoderm 3",
                                    "6"  = "Lateral plate mesoderm 1",
                                    "7"  = "Caudal epiblast/Early neural tube", 
                                    "8"  = "Dermomyotome", 
                                    "9"  = "Lateral plate mesoderm 2",
                                    "10" = "Extraembryonic mesoderm 1",
                                    "11" = "Surface ectoderm 3", 
                                    "12" = "Extraembryonic mesoderm 2",
                                    "13" = "Somitic mesoderm",
                                    "14" = "Definitive endoderm",
                                    "15" = "Hemato-endothelial progenitors 1",
                                    "16" = "Extraembryonic endoderm 2",
                                    "17" = "Surface ectoderm 2",
                                    "18" = "Amnion",
                                    "19" = "Neural crest",
                                    "20" = "Extraembryonic endoderm 1",
                                    "21" = "Mesenchyme 2",
                                    "22" = "Anterior primitive streak/Node",
                                    "23" = "Neurons",
                                    "24" = "Trophoblast",
                                    "25" = "Surface ectoderm 1",
                                    "26" = "Muscle",
                                    "27" = "Hemato-endothelial progenitors 2",
                                    "28" = "Unknown",
                                    "29" = "Early hematopoietic cells",
                                    "30" = "Extraembryonic endoderm 3",
                                    "31" = "PGCs",
                                    "32" = "Neural tube/Radial glia")
datasets.integrated$annotations.new <- Idents(datasets.integrated)
unique(datasets.integrated$annotations.new)

custom_order <- c(
  "Extraembryonic endoderm 1", "Extraembryonic endoderm 2", "Extraembryonic endoderm 3",
  "Amnion", "Trophoblast",
  "Extraembryonic mesoderm 1", "Extraembryonic mesoderm 2",
  "Epiblast/Neural plate", "Anterior primitive streak/Node", "Caudal epiblast/Early neural tube",
  "Nascent mesoderm/Presomitic mesoderm", "Somitic mesoderm", "Dermomyotome", "Muscle",
  "Lateral plate mesoderm 1", "Lateral plate mesoderm 2", "Lateral plate mesoderm 3",
  "Mesenchyme 1", "Mesenchyme 2", "Mesenchyme 3",
  "Hemato-endothelial progenitors 1", "Hemato-endothelial progenitors 2", "Early hematopoietic cells",
  "Definitive endoderm",
  "Neural tube/Radial glia", "Neurons", "Neural crest",
  "Surface ectoderm 1", "Surface ectoderm 2", "Surface ectoderm 3",
  "PGCs", "Unknown"
)

datasets.integrated$annotations.new.ordered <- factor(
  as.character(datasets.integrated$annotations.new),
  levels = custom_order
)

DimPlot(datasets.integrated, reduction = "reduction.new.params", group.by = "annotations.new.ordered")
saveRDS(datasets.integrated, "datasets.integrated.new.annotated.rds")

# Plot how well they agree with gast cell type predictions
meta <- datasets.integrated@meta.data
old_col <- "predictions.gast"         
new_col <- "annotations.new" 
ct <- table(meta[[old_col]], meta[[new_col]])

col_prop <- prop.table(ct, 2)
pheatmap(col_prop,
         color = colorRampPalette(c("white", "#E48C2A"))(100),
         cluster_rows = TRUE, cluster_cols = TRUE,
         display_numbers = TRUE, number_format = "%.2f")

