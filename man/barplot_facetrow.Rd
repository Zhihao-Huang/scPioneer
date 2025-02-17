% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tissue_similar.R
\name{barplot_facetrow}
\alias{barplot_facetrow}
\title{Barplot faceted in row. Plot data is from stat_meta_freq().}
\usage{
barplot_facetrow(
  df,
  color_palette = NULL,
  test.method = c("t.test", "wilcox.test"),
  adjust.method = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr",
    "none"),
  p.label = c("p.adj.signif", "p.adj", "p"),
  add.signif = T,
  return.stat = F
)
}
\description{
Barplot faceted in row. Plot data is from stat_meta_freq().
}
