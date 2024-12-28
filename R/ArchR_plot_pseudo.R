#' Merge pseudo paths on one plot
#' 
#' @export
ArchR_merge_pseudo <- function(object, traj_namelist, point.name,
                               color.line = NULL, color.point = NULL) {
  linelist <- lapply(traj_namelist, function(x) {
    p <- plotTrajectory(object, trajectory = x,
                        colorBy = "cellColData", name = x)
    dfline <- p[[1]]$layers[[3]]$data
    dfline$g <- x
    dfline
  })
  names(linelist) <- traj_namelist
  pathdf <- do.call(rbind, linelist)
  
  pointlist <- lapply(linelist, function(df) {
    df[c(1,nrow(df)),]
  })
  names(pointlist) <- traj_namelist
  dfp <- do.call(rbind, pointlist)
  dfp$g <- point.name
 
  df <- getEmbedding(object, embedding = 'UMAP', returnDF = TRUE)
  colnames(df) <- c('x','y')
  p <- ggplot(df, aes(x, y)) + geom_point(color = 'lightgrey', size = 0.5) 
  p <- p + geom_point(data = dfp, aes(x,y,fill = g),
                      shape = 21, color = 'black', size = 3)
  p <- p + geom_path(data = pathdf, aes(x, y, color = g),size = 0.8, 
                     arrow = arrow(type = "open", length = unit(0.1,"inches")))
  
  p <- p + ggrepel::geom_text_repel(data = dfp, aes(x,y, label = g),
                                    color = 'black', size = 4)
  p <- p + theme_bw()
  p <- p + theme(legend.position = 'none', 
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank())
  p <- p + xlab('UMAP_1') + ylab('UMAP_2')
  if (is.null(color.line)) {
    color.line <- scPalette2(length(traj_namelist))
  }
  if (is.null(color.point)) {
    color.point <- scPalette2(length(traj_namelist)*2)
  }
  p <- p + scale_color_manual(values = color.line)
  p <- p + scale_fill_manual(values = color.point)
  return(p)
}