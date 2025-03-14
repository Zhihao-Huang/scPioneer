% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_heatmap.R
\name{scale_color_heatmap}
\alias{scale_color_heatmap}
\title{Scale color of heatmap for pheatmap.}
\usage{
scale_color_heatmap(
  OR,
  colors = NULL,
  mid.point.positive = 1,
  mid.point.negative = 0,
  length.node.num = NULL,
  legend.breaks = NULL,
  cut.off = NULL,
  legend.title = NULL,
  paletteLength = 100
)
}
\arguments{
\item{OR}{matrix of value. The format is features x samples}

\item{mid.point.positive}{Value to indicate the middle point of legend. 
And all the value in matrix are positive (>0).}

\item{mid.point.negative}{Value to indicate the middle point of legend. 
And the minimum value in matrix is negative (<0).}
}
\description{
Scale color of heatmap for pheatmap.
}
