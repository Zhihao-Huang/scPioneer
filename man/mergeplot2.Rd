% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_merge.R
\name{mergeplot2}
\alias{mergeplot2}
\title{3 color plot of two co-expression genes.}
\usage{
mergeplot2(
  SeuratS4,
  marker_to_plot,
  active.assay = "RNA",
  reduction = "umap",
  expr.cut.off = NULL
)
}
\arguments{
\item{SeuratS4}{a Seurat object.}

\item{marker_to_plot}{two Genes to plot.}

\item{expr.cut.off}{Expression below Cut off will be assigned to 0.}
}
\description{
3 color plot of two co-expression genes.
}
\examples{
data('pbmc')
marker_to_plot <- c("CD4","LEF1")
mergeplot2(pbmc,marker_to_plot, expr.cut.off = 0.1)
}
