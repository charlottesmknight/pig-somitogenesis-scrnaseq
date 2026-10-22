# Needs atlas.subset.final.v2.rds (from script 6)
# Figures produced
#   - UMAP by cell type (annotations.v2), by dataset, and by embryonic day
#   - Feature plots of mesoderm markers (TBXT, TBX6, PAX3, MEOX1, DMRT2, MYOD1)
#   - Signalling module score dotplot
#   - Marker module heatmap, all clusters (Suppl. Fig. 1)
#   - Ridgeline of marker gene expression across embryonic day
#   - Ridgeline of population proportions across embryonic day

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(grid)
library(ggridges)

# V2 annotations with cardiac muscle cluster separated
asf <- readRDS("atlas.subset.final.v2.rds")
DefaultAssay(asf) <- "RNA"
asf[["RNA"]] <- JoinLayers(asf[["RNA"]])
Idents(asf) <- "annotations.v2"
rna_genes <- rownames(asf[["RNA"]])

# Colours
subset_palette <- c(
  "Epiblast 1"                        = "#AFD1CEFF",
  "Epiblast 2"                        = "#7CC3C9FF",
  "Primitive streak"                  = "#75AAA6FF",
  "Caudal epiblast"                   = "#519B9AFF",
  "Neural tube/Radial glia"           = "#1D7069FF",
  "Neurons"                           = "#03424BFF",
  "Anterior primitive streak"         = "#AFC7D6",
  "Node/Notochord"                    = "#1a5270",
  "Definitive endoderm"               = "#508ca7",
  "Nascent mesoderm"                  = "#EFBD9C",
  "Presomitic mesoderm"               = "#E69158",
  "Somitic mesoderm"                  = "#E1650BFF",
  "Dermomyotome"                      = "#BB4604",
  "Myoblasts"                         = "#A83600",
  "Myocytes/Myotubes"                 = "#802900",
  "Myofibers"                         = "#601F00",
  "Sclerotome"                        = "#82533A",
  "Somitic connective mesenchyme"     = "#996244FF",
  "Head mesoderm"                     = "#D0AB95",
  "Nascent mesoderm/Cardiac mesoderm" = "#C05F1F"
)
subset_palette_v2 <- c(subset_palette, "Cardiac muscle" = "#7A3E65")

stage_cols <- c(
  "E11.5" = "#DDEEF5", "E12"   = "#C3E1ED", "E12.5" = "#A9D5E5",
  "E13"   = "#8FC8DD", "E13.5" = "#73AEC6", "E14"   = "#5894AE",
  "E14.5" = "#3E7795", "E15"   = "#285D7E", "E16"   = "#1C5372",
  "E18"   = "#154263", "E21"   = "#113457", "E28"   = "#0E2A4D"
)

dataset_cols <- c("#C3E1ED", "#5894AE")

heatmap_cols <- c(
  Epiblast                         = "#AFD1CEFF",
  Primitive_Streak_Caudal_Epiblast = "#75AAA6FF",
  Anterior_Primitive_Streak        = "#AFC7D6",
  Node_Notochord                   = "#1a5270",
  Definitive_Endoderm              = "#508ca7",
  Nascent_Mesoderm                 = "#EFBD9C",
  Lateral_Plate_Cardiac_Mesoderm   = "#E46C2AFF",
  Head_Mesoderm                    = "#D0AB95",
  Presomitic_Mesoderm              = "#E69158",
  Somitic_Mesoderm                 = "#E1650BFF",
  Dermomyotome                     = "#BB4604",
  Sclerotome                       = "#82533A",
  Connective_Mesenchyme            = "#996244FF",
  Myoblasts                        = "#A83600",
  Myocytes_Myotubes                = "#802900",
  Myofibers                        = "#601F00",
  Cardiac_Muscle                   = "#7A3E65",
  Neural_Tube                      = "#1D7069FF",
  Neurons                          = "#03424BFF"
)

hm_cols   <- colorRampPalette(c("#5BA3BC", "white", "#CD4F38FF"))(100)
hm_breaks <- seq(-2, 2, length.out = 101)

# Cluster orders
reduced_clusters_asf_2 <- c(
  "Epiblast 1", "Epiblast 2", "Primitive streak", "Caudal epiblast",
  "Anterior primitive streak", "Node/Notochord", "Definitive endoderm",
  "Nascent mesoderm", "Nascent mesoderm/Cardiac mesoderm", "Head mesoderm",
  "Presomitic mesoderm", "Somitic mesoderm", "Dermomyotome", "Sclerotome"
)

# UMAP legend order
custom_order_asf <- c(
  "Epiblast 1", "Epiblast 2", "Caudal epiblast", "Primitive streak",
  "Anterior primitive streak", "Node/Notochord", "Definitive endoderm",
  "Nascent mesoderm", "Nascent mesoderm/Cardiac mesoderm", "Head mesoderm",
  "Presomitic mesoderm", "Somitic mesoderm", "Dermomyotome", "Sclerotome",
  "Somitic connective mesenchyme", "Myoblasts", "Myocytes/Myotubes",
  "Myofibers", "Neural tube/Radial glia", "Neurons"
)

custom_order_v2 <- append(custom_order_asf, "Cardiac muscle",
                          after = which(custom_order_asf == "Myofibers"))

# marker heatmap column order
my_cluster_order <- c(
  "Epiblast 1", "Epiblast 2", "Primitive streak", "Caudal epiblast",
  "Anterior primitive streak", "Node/Notochord", "Definitive endoderm",
  "Nascent mesoderm", "Nascent mesoderm/Cardiac mesoderm", "Head mesoderm",
  "Presomitic mesoderm", "Somitic mesoderm", "Dermomyotome", "Sclerotome",
  "Somitic connective mesenchyme", "Myoblasts", "Myocytes/Myotubes",
  "Myofibers", "Neural tube/Radial glia", "Neurons"
)

my_cluster_order_v2 <- append(my_cluster_order, "Cardiac muscle",
                              after = which(my_cluster_order == "Myofibers"))

# Marker modules (heatmap)
modules_merged_reduced <- list(
  Epiblast = c("OTX2", "DNMT3B", "NANOG"),
  Primitive_Streak_Caudal_Epiblast = c("CDX2", "CDX1", "FGF8", "GBX2"),
  Anterior_Primitive_Streak = c("NODAL", "GSC", "MIXL1"),
  Node_Notochord = c("SHH", "CHRD"),
  Definitive_Endoderm = c("SOX17", "CER1"),
  Nascent_Mesoderm = c("MESP1", "TBXT", "WNT5A", "WNT3A"),
  Lateral_Plate_Cardiac_Mesoderm = c("KDR", "HAND1", "HAND2", "FOXF1"),
  Head_Mesoderm = c("ALX1", "TCF21"),
  Presomitic_Mesoderm = c("DLL1", "TBX6", "MSGN1", "HES7"),
  Somitic_Mesoderm = c("MEOX1", "RIPPLY1", "MESP2", "ALDH1A2"),
  Dermomyotome = c("DMRT2", "PAX3", "MEOX2", "SIX1"),
  Sclerotome = c("PAX1", "PAX9", "SOX9"),
  Connective_Mesenchyme = c("PRRX2", "PDGFRB"),
  Myoblasts = c("MYF5", "PAX7"),
  Myocytes_Myotubes = c("MYOG", "MYMK", "MYOD1"),
  Myofibers      = c("MYH3", "MYL1", "MYLPF", "TNNI2"),
  Cardiac_Muscle = c("NKX2-5", "MYBPC3", "NPPB"),
  Neural_Tube = c("PAX6", "HES5", "NKX6-1", "IRX3"),
  Neurons = c("NEUROD1", "DCX", "STMN2")
)

# Signalling modules (dotplot)
signalling_modules_final_2 <- list(
  WNT = c("AXIN2", "NKD1", "SP5", "NOTUM", "LEF1", "TCF7", "TCF7L2", "WNT3",
          "WNT3A", "WNT8A", "WNT5A", "ROR2", "LRP6", "RNF43", "ZNRF3",
          "APCDD1", "FZD5", "FZD7", "FZD8"),
  FGF_posterior = c("FGF8", "DUSP6", "SPRY2", "SPRY4", "ETV4", "ETV5", "FGFR1",
                    "FGFR2", "DUSP4", "DUSP5", "SPRY1", "SPRED1", "SPRED2",
                    "IL17RD", "ETV1", "FGF4", "FGF17", "CYP26A1", "FLRT3"),
  NODAL = c("NODAL", "GDF1", "GDF3", "TDGF1", "FOXH1", "PITX2", "LEFTY1",
            "LEFTY2", "ACVR1B", "SMAD2", "SMAD3", "SMAD7", "ACVR2A", "ACVR2B",
            "ACVR1C"),
  TGFB_ACTIVIN = c("SERPINE1", "SKIL", "INHBA", "INHBB", "TGFB1", "TGFB2",
                   "TGFB3", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD7",
                   "ACVR2A", "ACVR2B", "GDF11", "LTBP1"),
  BMP = c("ID1", "ID2", "ID3", "MSX1", "MSX2", "BMP2", "BMP4", "BMPR1A",
          "BMPR2", "SMAD1", "SMAD5", "SMAD6", "BAMBI", "ID4", "BMP5", "BMP6",
          "BMP7", "ACVR1", "BMPR1B", "SMAD9", "GDF6"),
  RA_anterior = c("ALDH1A2", "RDH10", "STRA6", "CRABP2", "DHRS3", "RARB",
                  "ALDH1A1", "ALDH1A3", "RARA", "RARG", "RXRA", "RXRG",
                  "CRABP1", "RBP1", "FGF18"),
  NOTCH = c("DLL1", "DLL3", "HES7", "HES1", "LFNG", "HEY1", "HEY2", "NRARP",
            "NOTCH1", "NOTCH2", "JAG1", "RBPJ", "HEYL", "DTX1", "MFNG", "RFNG",
            "MIB1"),
  YAP_TAZ = c("CCN1", "CCN2", "ANKRD1", "AMOTL2", "THBS1", "AXL", "CRIM1"),
  SHH = c("GLI1", "PTCH1", "HHIP", "SHH", "GLI2", "SMO", "PTCH2", "GAS1",
          "CDON", "BOC", "DISP1", "SCUBE2", "HHAT", "SUFU", "GPR161")
)

signalling_order <- c("WNT", "FGF_posterior", "NODAL", "TGFB_ACTIVIN", "BMP",
                      "RA_anterior", "NOTCH", "YAP_TAZ", "SHH")

signalling_labels <- c(
  WNT = "WNT", FGF_posterior = "FGF", NODAL = "NODAL",
  TGFB_ACTIVIN = "TGFβ / Activin", BMP = "BMP",
  RA_anterior = "Retinoic acid", NOTCH = "NOTCH",
  YAP_TAZ = "YAP / TAZ", SHH = "SHH"
)

# Functions
umap_plot <- function(group, cols, ncol = NULL) {
  DimPlot(asf, reduction = "umap_test", group.by = group, cols = cols,
          shuffle = TRUE, raster = FALSE, pt.size = 1) +
    coord_fixed() +
    labs(x = "UMAP 1", y = "UMAP 2") +
    guides(color = guide_legend(
      ncol = ncol,
      override.aes = list(shape = 15, size = 6, stroke = 0))) +
    theme(plot.title      = element_blank(),
          legend.title    = element_text(size = 14, face = "bold"),
          legend.text     = element_text(size = 14),
          legend.key.size = unit(0.6, "cm"))
}

marker_heatmap <- function(modules, cluster_order, avg_fun, group_by = "annotations.v2",
                           col_labels = NULL) {
  mods  <- lapply(modules, intersect, y = rna_genes)
  mods  <- mods[lengths(mods) > 0]
  genes <- unique(unlist(mods))
  
  avg   <- avg_fun(asf, assays = "RNA", features = genes,
                   group.by = group_by)[["RNA"]]
  genes <- genes[genes %in% rownames(avg)]
  avg   <- avg[genes, cluster_order[cluster_order %in% colnames(avg)], drop = FALSE]
  
  if (!is.null(col_labels)) colnames(avg) <- col_labels[colnames(avg)]
  
  g2m      <- stack(mods)
  g2m      <- g2m[!duplicated(g2m$values), ]
  row_anno <- data.frame(Module = g2m$ind, row.names = g2m$values)
  row_anno <- row_anno[rownames(avg), , drop = FALSE]
  mod_vec  <- as.character(row_anno$Module)
  
  pheatmap(avg, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE,
           annotation_row = row_anno,
           annotation_colors = list(Module = heatmap_cols),
           gaps_row = which(mod_vec[-length(mod_vec)] != mod_vec[-1]),
           annotation_legend = TRUE, angle_col = 45,
           fontsize_row = 10, fontsize_col = 14, border_color = NA,
           color = hm_cols, breaks = hm_breaks, silent = TRUE)
}

score_modules <- function(obj, modules, prefix) {
  mods <- lapply(modules, function(g) unique(g[g %in% rownames(obj)]))
  obj  <- AddModuleScore(obj, features = mods, name = paste0(prefix, "__"),
                         seed = 42)
  old  <- paste0(prefix, "__", seq_along(mods))
  new  <- paste0(prefix, "__", names(mods))
  colnames(obj@meta.data)[match(old, colnames(obj@meta.data))] <- new
  print(lengths(mods))
  obj
}

module_long <- function(obj, prefix, modules, lab, row_order) {
  cols <- paste0(prefix, "__", names(modules))
  obj@meta.data %>%
    mutate(cluster = factor(annotations.v1, levels = reduced_clusters_asf_2)) %>%
    select(cluster, all_of(cols)) %>%
    pivot_longer(-cluster, names_to = "module", values_to = "score") %>%
    mutate(module = sub(paste0("^", prefix, "__"), "", module)) %>%
    group_by(cluster, module) %>%
    summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
    group_by(module) %>%
    mutate(score_scaled = (mean_score - min(mean_score)) /
             (max(mean_score) - min(mean_score))) %>%
    ungroup() %>%
    mutate(module = factor(lab[module], levels = rev(lab[row_order])))
}

module_dotplot <- function(df, heights = c(15, 1)) {
  
  p <- ggplot(df, aes(x = cluster, y = module)) +
    geom_hline(yintercept = seq_along(levels(df$module)),
               colour = "grey92", linewidth = 0.25) +
    geom_vline(xintercept = seq_along(levels(df$cluster)),
               colour = "grey92", linewidth = 0.25) +
    geom_point(aes(size = score_scaled, colour = score_scaled)) +
    scale_colour_gradientn(
      colours = c("#F2EDE4", "#B5482E"),
      limits  = c(0, 1),
      breaks  = c(0, 0.5, 1),
      labels  = c("0", "0.5", "1"),
      name    = "Module score",
      guide   = guide_colourbar(
        barwidth = unit(0.4, "cm"), barheight = unit(3, "cm"),
        frame.colour = "grey30", frame.linewidth = 0.3,
        ticks.colour = "grey30", ticks.linewidth = 0.3, order = 1)) +
    scale_size_continuous(
      range  = c(0.8, 7),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      name   = NULL,
      guide  = guide_legend(
        override.aes = list(colour = "grey35"), order = 2,
        keywidth = unit(0.5, "cm"), keyheight = unit(0.5, "cm"))) +
    scale_x_discrete(expand = expansion(add = 0.7)) +
    scale_y_discrete(expand = expansion(add = 0.7)) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 11) +
    theme(
      panel.border     = element_rect(colour = "grey25", fill = NA, linewidth = 0.5),
      panel.grid       = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.text.y      = element_text(size = 10, colour = "grey15"),
      axis.title       = element_blank(),
      axis.ticks.y     = element_line(colour = "grey40", linewidth = 0.3),
      axis.ticks.length = unit(0.12, "cm"),
      legend.title     = element_text(size = 9, colour = "grey15", lineheight = 1.1),
      legend.text      = element_text(size = 8.5, colour = "grey25", hjust = 0),
      legend.key       = element_blank(),
      legend.spacing.y = unit(0.3, "cm"),
      legend.position  = "right",
      legend.box       = "vertical",
      legend.box.just  = "left",
      plot.margin      = margin(t = 6, r = 6, b = 2, l = 6))
  
  strip_df <- data.frame(cluster = factor(levels(df$cluster),
                                          levels = levels(df$cluster)))
  
  strip <- ggplot(strip_df, aes(x = cluster, y = 1, fill = cluster)) +
    geom_tile(height = 1) +
    scale_fill_manual(values = subset_palette, guide = "none") +
    scale_x_discrete(expand = expansion(add = 0.7)) +
    theme_void() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9.5,
                                 colour = "grey15", margin = margin(t = 4)),
      plot.margin = margin(t = 2, r = 6, b = 6, l = 6))
  
  p / strip + plot_layout(heights = heights)
}


# UMAP - cell type
asf$annotations.v2.ordered <- factor(as.character(asf$annotations.v2),
                                     levels = custom_order_v2)
setdiff(unique(as.character(asf$annotations.v2)), custom_order_v2)   # expect character(0)

p_umap_type <- umap_plot("annotations.v2.ordered", subset_palette_v2, ncol = 2)
p_umap_type

# UMAP - dataset
asf$dataset <- recode(asf$dataset,
                      "pig_gastrulation_dataset" = "Gastrulation atlas",
                      "pig_myogenesis_data"      = "Myogenesis dataset")

p_umap_dataset <- umap_plot("dataset", dataset_cols)
p_umap_dataset

# UMAP - stage
asf$stage.simple <- sub("^((E[0-9]+(\\.[0-9]+)?)).*$", "\\1", asf$orig.ident)
asf$stage.simple <- factor(asf$stage.simple, levels = names(stage_cols))

p_umap_stage <- umap_plot("stage.simple", stage_cols)
p_umap_stage

# Feature plots
meso <- c("TBXT", "TBX6", "PAX3", "MEOX1", "DMRT2", "MYOD1")

p_feat <- FeaturePlot(asf, features = meso, reduction = "umap_test",
                      order = TRUE, min.cutoff = "q50", pt.size = 0.5, ncol = 3) &
  scale_color_gradient(low = "grey90", high = "#B5482E") &
  coord_fixed() &
  labs(x = "UMAP 1", y = "UMAP 2") &
  theme_bw() &
  theme(aspect.ratio      = 1,
        panel.grid        = element_blank(),
        plot.title        = element_text(face = "italic", size = 12, hjust = 0.5),
        axis.title        = element_text(size = 10),
        axis.text         = element_text(size = 8),
        legend.key.width  = unit(0.25, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.title      = element_text(size = 8),
        legend.text       = element_text(size = 7))
p_feat

# Signalling dotplot
asf_red <- subset(asf, idents = reduced_clusters_asf_2)
asf_red <- score_modules(asf_red, signalling_modules_final_2, "sig")
df_sig <- module_long(asf_red, "sig", signalling_modules_final_2,
                      signalling_labels, signalling_order)
p_sig <- module_dotplot(df_sig, heights = c(20, 1))
p_sig

# Marker heatmap - all clusters (Suppl. Fig. 1)
short_labels <- c("Epi 1", "Epi 2", "PS", "CE", "APS", "Node/Noto", "DE", "NM",
                  "NM/Cardiac", "Head Meso", "PSM", "Somitic Meso", "Dermo",
                  "Sclero", "SCM", "Myob", "Myoc/Myot", "Myof", "Cardiac",
                  "Neural T/RG", "Neur")
names(short_labels) <- my_cluster_order_v2

ph_markers <- marker_heatmap(modules_merged_reduced, my_cluster_order_v2, AggregateExpression,
                             col_labels = short_labels)
grid.newpage(); grid.draw(ph_markers$gtable)

# Ridgeline - marker expression across embryonic day (E11.5 to E21 only)
gene_sets <- list(
  posterior  = c("TBXT", "NKX1-2", "WNT3A"),
  psm        = c("MSGN1", "TBX6", "HES7"),
  boundary   = c("MESP2", "RIPPLY1"),
  somite     = c("TCF15", "PAX3", "MEOX1", "UNCX", "TBX18"),
  dermo      = c("MEOX2", "DMRT2", "PAX7"),
  sclerotome = c("PAX1", "PAX9"),
  myogenic   = c("MYF5", "MYOD1", "MYOG")
)

lapply(gene_sets, setdiff, y = rna_genes)      # genes missing from the object
gene_sets  <- lapply(gene_sets, intersect, y = rna_genes)
gene_order <- unlist(gene_sets, use.names = FALSE)
grp_of     <- setNames(rep(names(gene_sets), lengths(gene_sets)), gene_order)

grp_cols <- c(posterior  = subset_palette[["Nascent mesoderm"]],
              psm        = subset_palette[["Presomitic mesoderm"]],
              boundary   = subset_palette[["Presomitic mesoderm"]],
              somite     = subset_palette[["Somitic mesoderm"]],
              dermo      = subset_palette[["Dermomyotome"]],
              sclerotome = subset_palette[["Sclerotome"]],
              myogenic   = subset_palette[["Myoblasts"]])

stage_num <- data.frame(
  stage = as.numeric(sub("^E", "", as.character(asf$stage.simple)))
)

expr <- FetchData(asf, vars = gene_order, layer = "data")   # log-normalised

pct_stage <- data.frame(expr, stage_num, check.names = FALSE) %>%
  filter(stage <= 21) %>%                                                        # drop E28 before scaling
  pivot_longer(all_of(gene_order), names_to = "gene", values_to = "expr") %>%
  group_by(gene, stage) %>% summarise(avg = mean(expr), .groups = "drop") %>%
  group_by(gene) %>% mutate(rel = avg / pmax(max(avg), 1e-9)) %>% ungroup() %>%  # each gene scaled to its own peak
  mutate(gene = factor(gene, levels = gene_order),
         grp  = factor(grp_of[as.character(gene)], levels = names(gene_sets)))

stage_vals <- sort(unique(pct_stage$stage))

p_ridge_markers <- ggplot(pct_stage, aes(x = stage, y = gene, height = rel, fill = grp)) +
  geom_ridgeline(scale = 0.9, colour = "grey20", linewidth = 0.3, alpha = 0.9) +
  geom_vline(xintercept = 15.5, linetype = "dashed", linewidth = 0.3) +          # dataset boundary (E15/E16)
  scale_fill_manual(values = grp_cols, guide = "none") +
  scale_y_discrete(limits = rev(gene_order)) +
  scale_x_continuous(breaks = stage_vals, labels = paste0("E", stage_vals),
                     expand = c(0, 0)) +
  labs(x = "Embryonic day", y = NULL) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_text(face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 1))
p_ridge_markers

# Ridgeline - population proportions across embryonic day (E11.5 to E21 only)
ridge_pops <- c("Nascent mesoderm", "Nascent mesoderm/Cardiac mesoderm", "Head mesoderm",
                "Presomitic mesoderm", "Somitic mesoderm", "Dermomyotome", "Sclerotome",
                "Myoblasts", "Myocytes/Myotubes", "Myofibers",
                "Neural tube/Radial glia", "Neurons")

smp_pop <- data.frame(
  sample     = asf$orig.ident,
  stage      = as.numeric(sub("^E", "", as.character(asf$stage.simple))),
  population = as.character(asf$annotations.v2)
)

# proportion of each population per sample, averaged per stage
pop_stage <- smp_pop %>%
  filter(stage <= 21) %>%                                                   
  count(sample, stage, population) %>%
  complete(nesting(sample, stage), population, fill = list(n = 0)) %>%
  group_by(sample) %>% mutate(prop = n / sum(n)) %>% ungroup() %>%
  filter(population %in% ridge_pops) %>%
  group_by(population, stage) %>% summarise(mean_prop = mean(prop), .groups = "drop") %>%
  group_by(population) %>% mutate(rel = mean_prop / pmax(max(mean_prop), 1e-9)) %>% ungroup() %>%
  mutate(population = factor(population, levels = ridge_pops))

pop_stage_vals <- sort(unique(pop_stage$stage))

p_ridge_pops <- ggplot(pop_stage, aes(x = stage, y = population, height = rel, fill = population)) +
  geom_ridgeline(scale = 0.9, colour = "grey20", linewidth = 0.3, alpha = 0.9) +
  geom_vline(xintercept = 15.5, linetype = "dashed", linewidth = 0.3) +
  scale_fill_manual(values = subset_palette_v2, guide = "none") +
  scale_y_discrete(limits = rev(ridge_pops)) +
  scale_x_continuous(breaks = pop_stage_vals,
                     labels = paste0("E", pop_stage_vals),
                     expand = c(0, 0)) +
  labs(x = "Embryonic day", y = NULL) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
p_ridge_pops

