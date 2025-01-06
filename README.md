R package for preprocessing and visualization of scRNA-seq data friendly.

Install required packages:

install.packages('remotes')

remotes::install_version('Seurat', version = '5.0.1')

remotes::install_github('satijalab/seurat-wrappers')

install.packages("https://cran.r-project.org/src/contrib/Archive/SeuratObject/SeuratObject_5.0.1.tar.gz", repos=NULL, type="source")

deps <- c('ggplot2', 'Seurat', 'reshape2', 'dplyr', 'BiocParallel', 'clustree', 
          'RColorBrewer',  'pheatmap', 'SingleCellExperiment', 'gridExtra', 
          'network', 'GGally', 'png', 'org.Hs.eg.db', 
          'reticulate', 'progress', 'pbapply',  'circlize', 
          'cowplot', 'dendextend', 'ggsci',  
          'openxlsx', 'S4Vectors', 
          'gridBase', 'BiocNeighbors', 'limma', 'scDblFinder')

pkgs <- installed.packages(.libPaths())

for (i in deps) if (!i %in% rownames(pkgs)) BiocManager::install(i,update = F)

Then install scPioneer:

remotes::install_github('Zhihao-Huang/scPioneer')
