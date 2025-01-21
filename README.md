# Introduction
>ScPioneer is a R package for preprocessing and visualization of scRNA-seq data friendly.

# Installation

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


# Quick-start

PHASE 1: QC and clustering
```
### Create a samplelist including samplename and datadir.
df <- data.frame(samplename = 'PBMC', datadir = './data-raw/filtered_gene_bc_matrices/hg19/')
write.table(df, file = 'samplelist.txt', sep = '\t')

### load data
library(scPioneer)

### set arguments
param <- PHASE1_run_Seurat_v5_QC_clustering_param_template()
param$samplelist <- './samplelist.txt'
param$outdir <- './result/'
param$min.features <- 200
param$max.mt <- 10
param$min.cells <- 3
param$species <- 'Human'
param$normalize_method <- 'SCT'
param$detect.doublet <- "scDblFinder"
param$return.plot <- T

### run
resultlist <- PHASE1_run_Seurat_v5_QC_clustering(param)

names(resultlist)
# object plotlist
patchwork::wrap_plots(resultlist$plotlist)
```

<p align="center">
  <img width="750"  src="https://github.com/Zhihao-Huang/scPioneer/blob/main/data-raw/phage1.png">
</p>

PHASE 2: cell-type annotation
```
### Load preprecessed object from PHASE 1.
pbmc <- resultlist$object

### Perform annotation by SingleR
obj <- annocell(pbmc, species = 'Human', method = 'SingleR', raw_cluster = 'seurat_clusters')
p1 <- DimPlot_idx(obj) + ggtitle('singleR')

### Perform annotation by top markers
markerdf <- data.frame(celltypes = c('T','T','NK','NK','Mono', 'Mono','DC','DC','B','B','Platelet'), 
 markers = c('CD3D','CD3E','NCAM1','NKG7','CD14','FCGR3A','CST3','CD1C','CD79A','MS4A1','PPBP'))
colnames(markerdf)
# "celltypes" "markers"
obj <- annocell(pbmc, species = 'Human', method = 'topgene', markerdf = markerdf, raw_cluster = 'seurat_clusters')
p2 <- DimPlot_idx(obj) + ggtitle('top gene')

### Perform annotation by LLM model (OpenAI)
Idents(pbmc) <- pbmc$seurat_clusters
markers <- FindAllMarkers(pbmc, logfc.threshold = 0.5, test.use = 'MAST', only.pos = T)
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
obj <- annocell(pbmc, species = 'Human', method = 'angrycell', db = 'openai',
          DE = top10, raw_cluster = 'seurat_clusters',model = "gpt-3.5-turbo", seed = 1234,
          base_url = ,
          api_key = )
p3 <- DimPlot_idx(obj) + ggtitle('gpt-3.5-turbo')

### Perform annotation by LLM model (ollama. Local model, less accurate than OpenAI)
annodf <- anno_allama(top10)

patchwork::wrap_plots(resultlist$plotlist)
```

<p align="center">
  <img width="750"  src="https://github.com/Zhihao-Huang/scPioneer/blob/main/data-raw/phage2.png">
</p>

