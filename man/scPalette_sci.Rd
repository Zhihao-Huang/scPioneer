% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/color.R
\name{scPalette_sci}
\alias{scPalette_sci}
\title{Generate colors inspired by R package: ggsci.}
\usage{
scPalette_sci(
  n,
  palette = c("npg", "aaas", "nejm", "lancet", "jama", "jco", "ucscgb", "d3",
    "locuszoom", "igv", "uchicago", "startrek", "tron", "futurama", "rickandmorty",
    "simpsons", "gsea", "material"),
  ...
)
}
\arguments{
\item{n}{Number of colors. Max number of colors per palette is 10.}
}
\value{
A color palette for plotting
}
\description{
Generate colors inspired by R package: ggsci.
}
\examples{
scPalette_sci(10, palette = 'npg')

}
