% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_violin.R
\name{vioplot}
\alias{vioplot}
\title{scater plot of gene expression.}
\usage{
vioplot(
  data,
  paths,
  ident,
  select.celltypes = NULL,
  facet.toward = "row",
  do.log = F,
  strip.text.size = 8,
  strip.text.angle = NULL,
  strip.text.hjust = NULL,
  strip.text.vjust = NULL,
  panel.spacing = 0.5,
  clip = "on",
  color.use = NULL,
  jitter = F,
  cell.order = NULL,
  add.ave.point = T,
  ave.point.size = 1.5,
  add.line = T,
  line.size = 0.5,
  add.box = F,
  strip.text.switch = NULL,
  border.color = "black",
  border.size = 0.5,
  fill.expression = F,
  fill.color = c("#FEE090", "#D73027", "#A50026"),
  mid.color.point = "auto",
  legend.position = "top",
  legend.name = "",
  jitter.size = 0.1,
  text.x.size = 8,
  text.y.size = 8,
  hjust = 0,
  vjust = 0,
  legend.scale = NULL,
  legend.title.size = 8,
  legend.text.size = 8,
  legend.key.size = 1,
  axis.text.x.angle = NULL,
  remove.axis.x = F,
  remove.axis.y = F
)
}
\arguments{
\item{data}{a normalized expression matrix(gene x cell).}

\item{ident}{Cell annotation. Idents(seurat_object)}

\item{facet.toward}{'col' or 'row'.}

\item{strip.text.size}{Gene name text size.}

\item{gene}{Genes to plot.}
}
\value{
A ggplot object.
}
\description{
scater plot of gene expression.
}
\examples{
data('pbmc')
marker_to_plot <- c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A")
exp <- as.data.frame(t(as.matrix(pbmc@assays$RNA@data[marker_to_plot,, drop=FALSE])))
vioplot(exp, marker_to_plot, Idents(pbmc))
vioplot(exp, marker_to_plot, Idents(pbmc),fill.expression = T,add.ave.point = T,add.line = T,
jitter = T,add.box = T,legend.name = 'Average expression')
}
