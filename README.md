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
'''
### Automatically perform QC and clustering from raw rds.
param <- PHASE1_run_Seurat_v5_QC_clustering_param_template()
param$object <- 'PATH TO YOUR RAW RDS'
param$is_multidata <- 'TRUE'
param$sample_colname <- 'COLNAME OF YOUR SAMPLES'
param$detect.doublet <- 'scDblFinder'
param$outdir <- 'OUTPUT PATH'
obj <- PHASE1_run_Seurat_v5_QC_clustering(param)
'''
