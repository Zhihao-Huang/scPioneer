% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_highlight_cluster.R
\name{highallcl}
\alias{highallcl}
\title{highlight points of all cluster}
\usage{
highallcl(
  SeuratS4 = SeuratS4,
  reductions = "tsne",
  ncol = NULL,
  color = c("red", "grey")
)
}
\arguments{
\item{SeuratS4}{A Seurat object.}

\item{reduction}{Celect the reduction method.}
}
\value{
A ggplot objects.
}
\description{
highlight points of all cluster
}
\examples{
data('pbmc')
highallcl(pbmc,reduction = "umap")
}
