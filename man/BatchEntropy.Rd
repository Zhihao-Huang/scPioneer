% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BatchEntropy.R
\name{BatchEntropy}
\alias{BatchEntropy}
\title{Calculate the cell entropy according to the cell group. (From Liziyi, onco-fetal_utils.R)}
\usage{
BatchEntropy(
  input.data,
  group.id,
  k.used = 30,
  dimension.used = "raw",
  n.cores = 1
)
}
\arguments{
\item{input.data}{A matrix with each cell in a row. The column could be
tSNE coordinate, expression matrix or even a dist format.}

\item{group.id}{A vector with the same length as the row of input.data, or
a list of vectors.}

\item{k.used}{The k used to build a kNN-graph.}

\item{dimension.used}{The method to reduce the dimension, tSNE by default,
could also be PCA or raw.}
}
\value{
A vector with each cell's entropy.
}
\description{
Calculate the cell entropy according to the cell group. (From Liziyi, onco-fetal_utils.R)
}
\author{
Liziyi
}
