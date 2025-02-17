% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_heatmap_TCR_share.R
\name{shareTCR_heatmap}
\alias{shareTCR_heatmap}
\title{Heatmap of share TCR}
\usage{
shareTCR_heatmap(
  meta,
  group,
  annotation,
  clonotype_index,
  typemat,
  cell.order = NULL,
  keep.group = NULL,
  show.unique.TCR.group = NULL,
  legend.ncol = NULL,
  widths = c(1, 8, 6),
  x.size = 15,
  legend.heatmap.size = 10,
  legend.anno.size = 10,
  legend.heatmap.title.size = 15,
  legend.anno.title.size = 15,
  return.data = F
)
}
\description{
Heatmap of share TCR
}
\examples{
T_hgsoc <- meta200
typemat1 <- data.frame(Tissue = c('PT','PT'),
celltype = c('T08','T10'))
list1 <- shareTCR_heatmap(meta = T_hgsoc, group = 'Group_abb', annotation = 'Annotation_abb',
clonotype_index = 'clonotype_index', typemat = typemat1, return.data = T)
list1$plot

}
