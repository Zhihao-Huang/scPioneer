% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orthologs_mm_hs.R
\name{convert_symbol_mouse_human}
\alias{convert_symbol_mouse_human}
\title{Convert symbol gene names between mouse and human.}
\usage{
convert_symbol_mouse_human(
  geneset,
  species,
  version = 101,
  using.local.file = F,
  file.dir = "/BGFS1/home/huangzh/workspace/database/orthologs/mouse_human/",
  using.local.data = F,
  remove.duplicated = T,
  only.one.to.one = F
)
}
\arguments{
\item{geneset}{Vector of genes.}

\item{species}{Only human or mouse.}

\item{version}{Version of ensembl.}

\item{using.local.file}{Whether to use annotation file from local directory instead of connecting to internet.}

\item{file.dir}{Directory for annotation file.}

\item{remove.duplicated}{Whether to remove duplicated gene symbol.}

\item{only.one.to.one}{Whether to keep rows in which both gene symbol and ensembl ID are unique.}
}
\description{
Convert symbol gene names between mouse and human.
}
\examples{
df <- Get_orthologs_mouse_human(using.local.file=T)
filedir <- '/BGFS1/home/huangzh/workspace/database/orthologs/mouse_human/'
df <- Get_orthologs_mouse_human(version = 104, using.local.file = T,file.dir = filedir, only.one.to.one = T)

}
