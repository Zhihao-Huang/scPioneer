% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_gene_score.R
\name{boxplot_score_facetgene}
\alias{boxplot_score_facetgene}
\title{Box plot of gene score for each gene.}
\usage{
boxplot_score_facetgene(
  object,
  assay = "RNA",
  features,
  group.by = NULL,
  select.group = NULL,
  color.use = NULL,
  test.use = c("none", "t", "wilcox"),
  adjust.method = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
    "none"),
  plabel = c("p.adj.signif", "p.adj", "p"),
  ref.group = NULL,
  remove.bracket = F,
  adjust_p_text_gap = 1,
  adjust_p_text_height = 1,
  ylim = NULL,
  ylim.adjust = 1.1,
  return.stat = F,
  ...
)
}
\description{
Box plot of gene score for each gene.
}
\examples{
cellt <- c('Naive CD4 T', 'Memory CD4 T','CD8 T','NK')
subs <- subset(pbmc, Annotation \%in\% cellt)
geneset <- c('CCR7','TCF7','SELL')
boxplot_score_facetgene(subs, features = geneset)

}
