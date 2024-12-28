# AI friendly plot----
options(bitmapType = 'cairo')
deepcopy <- function(p) {
  unserialize(serialize(p, NULL))
}

gg0point <- function(plot){
  newplot <- deepcopy(plot)
  newplot$layers[[1]]$aes_params$size <- -1
  return(newplot)
}

#' Raster points in UMAP in order to make umap pdf easily be loaded by AI.
#' This function is suitable for plot with single coordinates.
#' 
#' @author Li Ziyi
#' @export 
ggAIplot <- function(plot, calibration = c(0,0,0,0), width = 10, height = 10, dpi = 300){
  plot.png <- plot + cowplot::theme_nothing()
  ggsave(plot.png, filename = "./temp.png", width = width, height = height, dpi = dpi)
  img <- png::readPNG("./temp.png")
  file.remove("./temp.png")
  blank.plot <- gg0point(plot)
  range.values <- c(
    ggplot_build(plot = blank.plot)$layout$panel_params[[1]]$x.range,
    ggplot_build(plot = blank.plot)$layout$panel_params[[1]]$y.range
  )
  plot.pdf <- blank.plot +
    annotation_raster(img, 
                      xmin = range.values[1] + calibration[1], 
                      xmax = range.values[2] + calibration[2],
                      ymin = range.values[3] + calibration[3], 
                      ymax = range.values[4] + calibration[4])
  return(plot.pdf)
}

annotation_raster.grid <- 
  function (raster, xmin, xmax, ymin, ymax, interpolate = FALSE, data) {
    raster <- grDevices::as.raster(raster)
    layer(data = data, mapping = NULL, stat = StatIdentity, 
          position = PositionIdentity, geom = GeomRasterAnn, inherit.aes = FALSE, 
          params = list(raster = raster, xmin = xmin, xmax = xmax, 
                        ymin = ymin, ymax = ymax, interpolate = interpolate))
  }
#' Raster points in UMAP in order to make umap pdf easily be loaded by AI. 
#' This function is suitable for multi-facet plots.
#' 
#' @author Li Ziyi
#' @export 
ggAIplot.grid <- function(plot, facet.by, calibration = c(0,0,0,0), width = 10, height = 10, dpi = 300){
  plot.data <- plot$data
  plot.pdf <- gg0point(plot)
  range.values <- c(
    ggplot_build(plot = plot.pdf)$layout$panel_params[[1]]$x.range,
    ggplot_build(plot = plot.pdf)$layout$panel_params[[1]]$y.range
  )
  for(variable in unique(plot.data[,facet.by])){
    temp.plot.png <- plot + cowplot::theme_nothing()
    temp.plot.png$data <- temp.plot.png$data[temp.plot.png$data[,facet.by] == variable,]
    ggsave(temp.plot.png, filename = "./temp.png", width = width, height = height, dpi = dpi)
    img <- png::readPNG("./temp.png")
    file.remove("./temp.png")
    ggimg <- annotation_raster.grid(img, 
                                    xmin = range.values[1] + calibration[1], 
                                    xmax = range.values[2] + calibration[2],
                                    ymin = range.values[3] + calibration[3], 
                                    ymax = range.values[4] + calibration[4],
                                    data = temp.plot.png$data[1,])
    plot.pdf <- plot.pdf + ggimg
  }
  return(plot.pdf)
}
