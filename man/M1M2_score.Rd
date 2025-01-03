% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_gene_score.R
\name{M1M2_score}
\alias{M1M2_score}
\title{Define M1/M2 phenotypes of Myeloid cells as mean expression of gene signatures, and plot data.
M1/M2 markers were from DOI: 10.1016/j.cell.2020.03.048 and DOI: 10.1016/j.cell.2018.05.060}
\usage{
M1M2_score(
  object,
  Macrophage.names,
  return.data = F,
  M1M2list = NULL,
  AddModuleScore = F,
  plot.type = "violin",
  correlation_centre = T
)
}
\arguments{
\item{object}{A Seurat object}

\item{Macrophage.names}{A vector of all the Macrophage sub-type names.}

\item{return.data}{Return statistic data.}

\item{plot.type}{The type of plot. Only 'violin', 'box' or 'correlation' is provided.}
}
\value{
A ggplot object of violin, box or heatmap.
}
\description{
Define M1/M2 phenotypes of Myeloid cells as mean expression of gene signatures, and plot data.
M1/M2 markers were from DOI: 10.1016/j.cell.2020.03.048 and DOI: 10.1016/j.cell.2018.05.060
}
