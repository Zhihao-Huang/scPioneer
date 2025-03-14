% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DEG.R
\name{DEG_group}
\alias{DEG_group}
\title{Pairwise DEG of multi-groups per cell type using Seurat::FindMarkers and return DEGlist.}
\usage{
DEG_group(
  object,
  groupname,
  test.use = c("wilcox", "t", "MAST", "DESeq2", "limma-trend", "limma-voom", "MAST_RE",
    "glmmTMB"),
  min.pct = 0.25,
  logfc.threshold = 0.2,
  batch.colname = NULL,
  min.cell.num = 3,
  g.cores = 1,
  save.path = NULL,
  annotation = NULL,
  ...
)
}
\description{
Pairwise DEG of multi-groups per cell type using Seurat::FindMarkers and return DEGlist.
}
\examples{
pbmc$batch <- c(rep('sample1',1000), rep('sample2',1000),rep('sample3',638))
subs4 <- subset(pbmc, Annotation \%in\% c('NK','DC','B'))
DEGlist <- DEG_group(subs4, groupname = 'Annotation',test.use = 'limma-trend')
DEGlist <- DEG_group(subs4, groupname = 'Annotation',test.use = 'limma-trend',batch.colname = 'batch')
}
