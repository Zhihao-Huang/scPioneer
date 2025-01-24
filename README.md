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
          'limma', 'clusterProfiler','glmGamPoi')
exist_pkgs <- installed.packages(.libPaths())
for (i in deps) if (!i %in% rownames(exist_pkgs)) BiocManager::install(i,update = F)

remotes::install_github("satijalab/seurat-wrappers", "seurat5")

```


Then install scPioneer:

```
remotes::install_github('Zhihao-Huang/scPioneer')
```


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


# PHASE 2: cell-type annotation
Load data:
```
### Load preprecessed object from PHASE 1.
pbmc <- resultlist$object
```
Perform annotation by reference (SingleR)
```
### Perform annotation by reference (SingleR)
obj <- annocell(pbmc, species = 'Human', method = 'SingleR', raw_cluster = 'seurat_clusters')
p1 <- DimPlot_idx(obj) + ggtitle('singleR')
p1
```
<p align="center">
  <img width="250"  src="https://github.com/Zhihao-Huang/scPioneer/blob/main/data-raw/anno_singler.png">
</p>

Perform annotation by LLM model:
```
### Get DEGs by FindAllMarkers
Idents(pbmc) <- pbmc$seurat_clusters
markers <- FindAllMarkers(pbmc, logfc.threshold = 0.5, test.use = 'MAST', only.pos = T)
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
```
OpenAI GPT-3.5-turbo (Getting your api key from: https://openai.com/):
```
obj <- annocell(pbmc, species = 'Human', method = 'angrycell', db = 'openai',
                DE = top10, raw_cluster = 'seurat_clusters',model = "gpt-3.5-turbo", seed = 1234,
                base_url = "http://chatapi.littlewheat.com/v1",
                api_key = )
p2 <- DimPlot_idx(obj) + ggtitle('gpt-3.5-turbo')
p2
```
<p align="center">
  <img width="250"  src="https://github.com/Zhihao-Huang/scPioneer/blob/main/data-raw/anno_gpt-3.5-turbo.png">
</p>

DeepSeek R1 (Getting your api key from: https://api-docs.deepseek.com/):
```
obj <- annocell(pbmc, species = 'Human', method = 'angrycell', db = 'openai',
                DE = top10, raw_cluster = 'seurat_clusters',model = "deepseek-chat", seed = 1234,
                base_url = 'https://api.deepseek.com/v1',
                api_key = )
p3 <- DimPlot_idx(obj) + ggtitle('deepseek-chat')
p3
```
<p align="center">
  <img width="250"  src="https://github.com/Zhihao-Huang/scPioneer/blob/main/data-raw/anno_gpt-3.5-turbo.png">
</p>

Perform annotation by top markers:
```
### Perform annotation by top markers
markerdf <- data.frame(celltypes = c('T','T','NK','NK','Mono', 'Mono','DC','DC','B','B','Platelet'), 
                       markers = c('CD3D','CD3E','NCAM1','NKG7','CD14','FCGR3A','CST3','CD1C','CD79A','MS4A1','PPBP'))
colnames(markerdf)
# "celltypes" "markers"
obj <- annocell(pbmc, species = 'Human', method = 'topgene', markerdf = markerdf, raw_cluster = 'seurat_clusters')
p4 <- DimPlot_idx(obj) + ggtitle('top gene')
p4
```
<p align="center">
  <img width="250"  src="https://github.com/Zhihao-Huang/scPioneer/blob/main/data-raw/anno_topgene.png">
</p>

# SessionInfo
<details>
<summary>SessionInfo</summary>
```
sessionInfo()
```
R version 4.4.2 (2024-10-31)
Platform: x86_64-conda-linux-gnu
Running under: Rocky Linux 8.10 (Green Obsidian)

Matrix products: default
BLAS/LAPACK: /home/jasper/.conda/envs/scPioneer/lib/libopenblasp-r0.3.28.so;  LAPACK version 3.12.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

time zone: Asia/Hong_Kong
tzcode source: system (glibc)

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods
[8] base

other attached packages:
 [1] scCustomize_3.0.1  scPioneer_2.1.0    dbplyr_2.3.4       dplyr_1.1.4
 [5] Matrix_1.7-1       reshape2_1.4.4     Seurat_5.0.0       SeuratObject_5.0.2
 [9] sp_2.1-4           ggplot2_3.5.1

loaded via a namespace (and not attached):
  [1] IRanges_2.40.1              randomcoloR_1.1.0.1
  [3] R.methodsS3_1.8.2           progress_1.2.3
  [5] goftest_1.2-3               Biostrings_2.74.1
  [7] HDF5Array_1.34.0            vctrs_0.6.5
  [9] spatstat.random_3.3-2       digest_0.6.37
 [11] png_0.1-8                   shape_1.4.6.1
 [13] registry_0.5-1              gypsum_1.2.0
 [15] ggrepel_0.9.6               deldir_2.0-4
 [17] parallelly_1.41.0           MASS_7.3-64
 [19] MAST_1.32.0                 httpuv_1.6.15
 [21] foreach_1.5.2               BiocGenerics_0.52.0
 [23] withr_3.0.2                 ggrastr_1.0.2
 [25] survival_3.8-3              memoise_2.0.1
 [27] ggbeeswarm_0.7.2            janitor_2.2.1
 [29] systemfonts_1.1.0           ragg_1.3.3
 [31] zoo_1.8-12                  GlobalOptions_0.1.2
 [33] V8_6.0.0                    pbapply_1.7-2
 [35] R.oo_1.27.0                 prettyunits_1.2.0
 [37] rematch2_2.1.2              KEGGREST_1.46.0
 [39] promises_1.3.2              httr_1.4.7
 [41] restfulr_0.0.15             globals_0.16.3
 [43] fitdistrplus_1.2-1          rhdf5filters_1.18.0
 [45] rstudioapi_0.17.1           rhdf5_2.50.2
 [47] UCSC.utils_1.2.0            miniUI_0.1.1.1
 [49] generics_0.1.3              curl_6.1.0
 [51] S4Vectors_0.44.0            zlibbioc_1.52.0
 [53] ScaledMatrix_1.14.0         polyclip_1.10-7
 [55] ca_0.71.1                   GenomeInfoDbData_1.2.13
 [57] ExperimentHub_2.14.0        SparseArray_1.6.0
 [59] xtable_1.8-4                stringr_1.5.1
 [61] doParallel_1.0.17           S4Arrays_1.6.0
 [63] BiocFileCache_2.14.0        hms_1.1.3
 [65] GenomicRanges_1.58.0        irlba_2.3.5.1
 [67] colorspace_2.1-1            filelock_1.0.3
 [69] ROCR_1.0-11                 isoband_0.2.7
 [71] reticulate_1.40.0           spatstat.data_3.1-4
 [73] snakecase_0.11.1            magrittr_2.0.3
 [75] lmtest_0.9-40               later_1.4.1
 [77] viridis_0.6.5               lattice_0.22-6
 [79] glmGamPoi_1.18.0            spatstat.geom_3.3-4
 [81] future.apply_1.11.3         genefilter_1.88.0
 [83] scattermore_1.2             XML_3.99-0.18
 [85] scuttle_1.16.0              cowplot_1.1.3
 [87] matrixStats_1.4.1           RcppAnnoy_0.0.22
 [89] pillar_1.10.0               nlme_3.1-166
 [91] iterators_1.0.14            compiler_4.4.2
 [93] beachmat_2.22.0             RSpectra_0.16-2
 [95] stringi_1.8.4               TSP_1.2-4
 [97] lubridate_1.9.4             tensor_1.5
 [99] SummarizedExperiment_1.36.0 GenomicAlignments_1.42.0
[101] plyr_1.8.9                  crayon_1.5.3
[103] abind_1.4-8                 BiocIO_1.16.0
[105] scater_1.34.0               locfit_1.5-9.10
[107] bit_4.5.0.1                 codetools_0.2-20
[109] textshaping_0.4.1           BiocSingular_1.22.0
[111] paletteer_1.6.0             alabaster.ranges_1.6.0
[113] GetoptLong_1.0.5            plotly_4.10.4
[115] mime_0.12                   splines_4.4.2
[117] circlize_0.4.16             Rcpp_1.0.13-1
[119] fastDummies_1.7.4           sparseMatrixStats_1.18.0
[121] utf8_1.2.4                  blob_1.2.4
[123] clue_0.3-66                 BiocVersion_3.20.0
[125] listenv_0.9.1               DelayedMatrixStats_1.28.0
[127] openxlsx_4.2.7.1            tibble_3.2.1
[129] statmod_1.5.0               pkgconfig_2.0.3
[131] pheatmap_1.0.12             tools_4.4.2
[133] cachem_1.1.0                RSQLite_2.3.9
[135] viridisLite_0.4.2           DBI_1.2.3
[137] celldex_1.16.0              scDblFinder_1.20.0
[139] fastmap_1.2.0               scales_1.3.0
[141] ica_1.0-3                   Rsamtools_2.22.0
[143] AnnotationHub_3.14.0        patchwork_1.3.0
[145] ggprism_1.0.5               BiocManager_1.30.25
[147] dotCall64_1.2               RANN_2.6.2
[149] alabaster.schemas_1.6.0     SingleR_2.8.0
[151] coro_1.1.0                  farver_2.1.2
[153] mgcv_1.9-1                  yaml_2.3.10
[155] MatrixGenerics_1.18.0       rtracklayer_1.66.0
[157] cli_3.6.3                   purrr_1.0.2
[159] stats4_4.4.2                leiden_0.4.3.1
[161] lifecycle_1.0.4             uwot_0.2.2
[163] Biobase_2.66.0              bluster_1.16.0
[165] BiocParallel_1.40.0         annotate_1.84.0
[167] timechange_0.3.0            gtable_0.3.6
[169] rjson_0.2.23                ggridges_0.5.6
[171] progressr_0.15.1            parallel_4.4.2
[173] limma_3.62.1                jsonlite_1.8.9
[175] edgeR_4.4.1                 RcppHNSW_0.6.0
[177] seriation_1.5.7             bitops_1.0-9
[179] bit64_4.5.2                 xgboost_1.7.8.1
[181] Rtsne_0.17                  alabaster.matrix_1.6.1
[183] spatstat.utils_3.1-1        BiocNeighbors_2.0.1
[185] zip_2.3.1                   metapod_1.14.0
[187] alabaster.se_1.6.0          dqrng_0.4.1
[189] spatstat.univar_3.1-1       R.utils_2.12.3
[191] lazyeval_0.2.2              alabaster.base_1.6.1
[193] shiny_1.10.0                htmltools_0.5.8.1
[195] sctransform_0.4.1           rappdirs_0.3.3
[197] glue_1.8.0                  spam_2.11-0
[199] httr2_1.0.7                 XVector_0.46.0
[201] RCurl_1.98-1.16             scran_1.34.0
[203] gridExtra_2.3               igraph_2.1.2
[205] R6_2.5.1                    sva_3.54.0
[207] tidyr_1.3.1                 DESeq2_1.46.0
[209] SingleCellExperiment_1.28.1 forcats_1.0.0
[211] labeling_0.4.3              cluster_2.1.8
[213] Rhdf5lib_1.28.0             GenomeInfoDb_1.42.1
[215] ellmer_0.1.0                DelayedArray_0.32.0
[217] tidyselect_1.2.1            vipor_0.4.7
[219] SEtools_1.20.0              AnnotationDbi_1.68.0
[221] future_1.34.0               S7_0.2.0
[223] rsvd_1.0.5                  munsell_0.5.1
[225] KernSmooth_2.23-26          data.table_1.16.4
[227] htmlwidgets_1.6.4           ComplexHeatmap_2.18.0
[229] RColorBrewer_1.1-3          rlang_1.1.4
[231] spatstat.sparse_3.1-0       spatstat.explore_3.3-3
[233] beeswarm_0.4.0              sechm_1.14.0
<details>

