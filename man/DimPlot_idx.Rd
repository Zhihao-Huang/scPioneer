% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DimPlot_idx.R
\name{DimPlot_idx}
\alias{DimPlot_idx}
\title{DimPlot with index prefix.}
\usage{
DimPlot_idx(
  object,
  group.by = NULL,
  reduction = "umap",
  label.idx = T,
  label.idx.size = 5,
  repel = F,
  idx.sep = " - ",
  idx.in.key = F,
  prefix.index = NULL,
  cover.outlier = F,
  alter.group.for.outlier = "seurat_clusters",
  rev.outlier.order = F,
  legend.ncol = 2,
  legend.byrow = F,
  legend.title = NULL,
  legend.position = "right",
  legend.key.size = 5,
  legend.text.size = 8,
  AI.friendly = F,
  return.obj = F,
  combine = F,
  cols = NULL,
  ...
)
}
\description{
DimPlot with index prefix.
}
\examples{
p1 <- DimPlot(pbmc)
p2 <- DimPlot_idx(pbmc)
p1 + p2

}
