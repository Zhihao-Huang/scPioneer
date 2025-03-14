% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tissue_similar.R
\name{draw_circle}
\alias{draw_circle}
\title{circle plot of network to show tissue origin.}
\usage{
draw_circle(
  net_pos,
  meta,
  group,
  query.group.name,
  annotation,
  group.order = NULL,
  cell.order = NULL,
  color.sector = NULL,
  color.link = NULL,
  rev_origin = F,
  frac.cutoff = 0,
  only.max.tissue = F,
  alpha.link = 0.4,
  labels.y = 1.5,
  legend.title = "",
  legend.title.position = "topleft",
  legend.direction = c("horizontal", "vertical")[1],
  legend.adjust.x = unit(4, "mm"),
  legend.adjust.y = unit(4, "mm"),
  legend.adjust = c("left", "bottom"),
  legend.nrow = NULL,
  legend.lwd = 2,
  legend.fontsize = 10,
  margin.par = c(1, 1, 1, 1),
  track.height = 0.1,
  return.df = F
)
}
\arguments{
\item{margin.par}{Margin of plot. Order: bottom, left, top, right.}
}
\description{
circle plot of network to show tissue origin.
}
