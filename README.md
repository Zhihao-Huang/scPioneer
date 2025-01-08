>ScPioneer is a R package for preprocessing and visualization of scRNA-seq data friendly.

Install required packages:
```
install.packages('remotes')
install.packages('BiocManager')

remotes::install_github("satijalab/seurat-data", "seurat5")

deps <- c('R.utils','ggplot2', 'reshape2', 'dplyr', 'BiocParallel', 'clustree', 'RColorBrewer',
          'pheatmap', 'SingleCellExperiment', 'gridExtra', 'network', 'GGally', 'png',
          'org.Hs.eg.db','org.Mm.eg.db', 'reticulate', 'progress', 'pbapply', 'circlize', 'cowplot', 
          'dendextend', 'ggsci','openxlsx', 'S4Vectors', 'gridBase', 'BiocNeighbors',
          'limma', 'clusterProfiler')
exist_pkgs <- installed.packages(.libPaths())
for (i in deps) if (!i %in% rownames(exist_pkgs)) BiocManager::install(i,update = F)

remotes::install_github("satijalab/seurat-wrappers", "seurat5")

```


Then install scPioneer:

`remotes::install_github('Zhihao-Huang/scPioneer')`


Quick-start:
Processing single sample.

```
### Automatically perform QC and clustering from raw rds.
library(scPioneer)
param <- PHASE1_run_Seurat_v5_QC_clustering_param_template()
param$object <- 'PATH TO YOUR RAW RDS'
param$outdir <- 'OUTPUT PATH'
obj <- PHASE1_run_Seurat_v5_QC_clustering(param)
```

Processing multi-samples from filtered matrix of Cellranger.
```
### Automatically perform QC and clustering from matrix files.
library(scPioneer)
samplelist <- data.frame(samplename = c('sample1','sample2'),
                              datadir = c('/XXX/sample1/outs/filtered_feature_bc_matrix/','/XXX/sample2/outs/filtered_feature_bc_matrix/'))
write.table(samplelist, file = samplelist_path, sep = '\t')
param <- PHASE1_run_Seurat_v5_QC_clustering_param_template()
param$samplelist <- samplelist_path
param$is_multidata <- 'TRUE'
param$sample_colname <- 'COLNAME OF YOUR SAMPLES'
param$detect.doublet <- 'scDblFinder'
param$outdir <- 'OUTPUT PATH'
param$species <- 'Mouse'
obj <- PHASE1_run_Seurat_v5_QC_clustering(param)
```
Processing multi-samples from Seruat object.
```
### Automatically perform QC and clustering from raw rds.
library(scPioneer)
param <- PHASE1_run_Seurat_v5_QC_clustering_param_template()
param$object <- 'PATH TO YOUR RAW RDS'
param$is_multidata <- 'TRUE'
param$sample_colname <- 'COLNAME OF YOUR SAMPLES'
param$detect.doublet <- 'scDblFinder'
param$outdir <- 'OUTPUT PATH'
param$species <- 'Mouse'
obj <- PHASE1_run_Seurat_v5_QC_clustering(param)
```
