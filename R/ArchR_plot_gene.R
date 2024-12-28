#' scatter plot of ArchR processed data
#' 
#' @export
plotgene_ArchR <- function(object, markerGenes, 
                           pt.size = 1,
                           AI.friendly = F,
                           AI.script = '/BGFS1/projectdata/jasper/project/ovarian_cancer/shell/functions/ggAIplot.R',
                           ncol = ceiling(length(markerGenes)^0.5),
                           plot.margin = unit(c(0, 0, 0, 0), "cm"),
                           size.title = 10){
  p <- plotEmbedding(
    ArchRProj = object,
    colorBy = "GeneScoreMatrix",
    name = markerGenes,
    embedding = "UMAP",
    size = pt.size, 
    imputeWeights = getImputeWeights(object)
  )
  if (ncol == 1) {
    return(p)
  }
  #Rearrange for grid plotting
  p2 <- lapply(names(p), function(x){
    plot <- p[[x]]  +
      theme_ArchR(baseSize = 6.5) +
      theme(
        plot.margin = plot.margin,
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5, size = size.title)
      ) + xlab('') + ylab('') + ggtitle(x) + guides(color = 'none', fill = 'none')
    if (AI.friendly) {
      source(AI.script)
      plot <- ggAIplot(plot)
    }
    plot
  })

  do.call(cowplot::plot_grid, c(list(ncol = ncol),p2))
}

#' violin plot of ArchR processed data
#' 
#' @export
vioplot_ArchR <- function(object, features, GroupBy = 'Clusters', 
                          colorBy = "GeneScoreMatrix",
                          ncol = 1,
                          addBoxPlot = T,plot.data = T, return.data = F, ...){ 
  margin = theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
  plist <- lapply(features, function(gene) {
    p <- plotGroups(
      ArchRProj = object, 
      groupBy = GroupBy, 
      colorBy = colorBy, 
      name = gene,
      plotAs = "violin",
      alpha = 0.4,
      addBoxPlot = addBoxPlot,
      ...
    )
    p + margin
  })
  names(plist) <- features
  plist2 <- lapply(features[1:length(features)-1], function(gene) {
    plist[[gene]] + theme(axis.text.x = element_blank(), 
                          axis.title.x = element_blank(), 
                          axis.ticks.x = element_blank())
  })
  names(plist2) <- features[-length(features)]
  plist2[[features[length(features)]]] <- plist[[features[length(features)]]]
  if (plot.data) {
    do.call(gridExtra::grid.arrange, c(plist2, ncol = ncol))
  }
  if (return.data) {
    return(plist2)
  }
}
