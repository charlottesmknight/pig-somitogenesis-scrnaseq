# pig-somitogenesis-scrnaseq

**Knight et al. (2026).** *Species-specific WNT response shapes somite morphogenesis in mammalian embryo models.*
Integration and annotation of pig embryonic single-cell RNA-seq datasets (GEO: GSE236766, GSE206914).
Run scripts 1–7 in order in a single R session

## Scripts

### 1. Integration.R
- **Inputs:** `Pig_atlas_Simpson_2024.RDS` (pig gastrulation atlas, Simpson et al., 2024) and `Myo_data_Cai_2023.RDS` (pig myogenesis dataset, Cai et al., 2023).
- **Outputs:** `datasets.integrated`, `pig_gastrulation_dataset`, `pig_myogenesis_data`.

### 2. Label_Transfer.R
- **Inputs:** `datasets.integrated`, `pig_gastrulation_dataset`, `pig_myogenesis_data` (from script 1); `Marioni_Atlas_pig_names_V5.RDS` (mouse gastrulation atlas, Pijuan-Sala et al., 2019, with mouse gene names converted to pig gene names; produced separately).
- Reclusters the integrated object and transfers cell type labels from the gastrulation atlas, the myogenesis dataset and the mouse atlas.
- **Outputs:** `datasets.integrated` with `predictions.gast`, `predictions.myo`, `predictions.mouse` and `predictions.merged`.

### 3. Clustering.R
- **Inputs:** `datasets.integrated` (from script 2).
- Final PCA, clustering (resolution 0.8) and UMAP of the integrated object.
- **Outputs:** `datasets.integrated` with `clusters.new.params` and the `reduction.new.params` UMAP.

### 4. Annotation_Full_Dataset.R
- **Inputs:** `datasets.integrated` (from script 3).
- Identifies cluster markers, annotates the 33 clusters and compares the annotation with the transferred gastrulation atlas labels.
- **Outputs:** `datasets.integrated` with `annotations.new` and `annotations.new.ordered`; `markers_res0.8`; `datasets.integrated.new.annotated.rds`.

### 5. Subset_Annotation.R
- **Inputs:** `datasets.integrated` (from script 4).
- Subsets nine neural and mesodermal populations, reclusters them, scores cell type gene modules and annotates the clusters.
- **Outputs:** `atlas.subset` with `annotations.v1`.

### 6. Myofiber_Cluster_Split.R
- **Inputs:** `atlas.subset` (from script 5).
- Separates cardiac muscle from skeletal myofibres within the myofibre cluster using cardiac and skeletal module scores.
- **Outputs:** `atlas.subset.final.v2.rds` (with `annotations.v2`).

### 7. Plots_Fig1.R
- **Inputs:** `atlas.subset.final.v2.rds` (from script 6).
- Generates the Figure 1 and Supplementary Figure 1 panels: UMAPs (cell type, dataset, embryonic day), marker feature plots, signalling module dotplot, marker heatmap and ridgeline plots.

## Data availability
- Pig gastrulation atlas (Simpson et al., 2024): GEO [GSE236766]
- Pig myogenesis dataset (Cai et al., 2023): GEO [GSE206914]
- Mouse gastrulation atlas (Pijuan-Sala et al., 2019): ArrayExpress E-MTAB-6967.
- Input objects (`Pig_atlas_Simpson_2024.RDS`, `Myo_data_Cai_2023.RDS`) and processed objects (`datasets.integrated.new.annotated.rds`, `atlas.subset.final.v2.rds`): Zenodo [DOI].

## Software
- R 4.4.1
- Seurat 5.4.0
- SeuratObject 5.3.0
- Matrix 1.7.4
- dplyr 1.1.4
- tidyr 1.3.2
- tibble 3.3.0
- ggplot2 4.0.1
- patchwork 1.3.2
- pheatmap 1.0.13
- ggridges 0.5.7

## Citation


## License

MIT
