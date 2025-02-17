% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_angrycell_result.R
\name{plot_pvalue}
\alias{plot_pvalue}
\title{Plot adjust p value of cell types.}
\usage{
plot_pvalue(objs4, reduction = "umap")
}
\arguments{
\item{candidate}{Genes to test.}

\item{N}{Number of cells.}

\item{db}{Marker genes database for usage.}
}
\value{
ggplot object.
}
\description{
Plot adjust p value of cell types.
}
\examples{
anno <- angrycell(pbmc)
pbmc@meta.data$Orig_Idents <- Idents(pbmc)
pbmc@meta.data$Cellname <- rownames(pbmc@meta.data)
pbmc@meta.data <- inner_join(pbmc@meta.data,anno)
rownames(pbmc@meta.data) <- pbmc@meta.data$Cellname
plot_pvalue(pbmc)
}
