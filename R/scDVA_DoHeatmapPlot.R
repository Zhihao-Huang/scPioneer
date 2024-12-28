#' Do the geneset heatmap plot.
#'
#'
#' @param plot.data The plot data created from getPlotData function.
#' @param genes The genes used in the plot.
#' @param group.by The cells should be ordered by or grouped by.
#' @param color.panel The color palette used when color by genes expression.
#' @param font.size The size of text, default by 15.
#' @return A heatmap plot.
#' 
#' @export
DoHeatmapPlot <- function(plot.data, genes, group.by, color.panel, 
                          only.heatmap = F,
                          coord.flip = F,
                          font.size = 15){
  if(length(genes) > 1){
    expression.matrix <- t(plot.data[, genes])
    genes <- genes[rowSums(expression.matrix) != 0]
    expression.matrix <- expression.matrix[rowSums(expression.matrix) != 0, ]
    groupid <- factor(plot.data[, group.by])
    
    if(length(unique(groupid)) > 1){
      # Boxplot in the left panel
      signature.mean.expression <- colMeans(expression.matrix[genes,])
      boxplot.data <- data.frame(Expression = signature.mean.expression, Group = groupid)
      p_boxplot <- 
        ggplot(boxplot.data, aes(x = Group, y = Expression)) + 
        geom_boxplot(color = "#1979B5", fill = "#1979B5", 
                     outlier.colour = "black", outlier.shape = 1) + 
        theme_bw() +
        theme(legend.position = "NULL",
              axis.text.y = element_text(size = font.size * 0.8),
              axis.text.x = element_text(size = font.size),
              axis.title = element_blank()) +
        scale_x_discrete(limits = rev(levels(as.factor(groupid)))) +
        coord_flip()
      dat <- ggplot_build(p_boxplot)$data[[1]]
      dat$xsp <- 1/2 * (dat$xmax + dat$xmin) - 1/4 * (dat$xmax - dat$xmin)
      dat$xep <- 1/2 * (dat$xmax + dat$xmin) + 1/4 * (dat$xmax - dat$xmin)
      p_boxplot <- 
        p_boxplot + 
        geom_segment(data = dat,  aes(x = xmin, xend = xmax, y = middle, yend = middle), colour = "#2BA147", size = 2) +
        geom_segment(data = dat,  aes(x = xsp, xend = xep, y = ymin, yend = ymin), colour = "black", size = 2) +
        geom_segment(data = dat,  aes(x = xsp, xend = xep, y = ymax, yend = ymax), colour = "black", size = 2)
      # Heatmap in the right panel
      if (color.panel %in% c("RdYlBu", "RdYlGn", "Spectral", "RdBu")) {
        myColorPalette <- colorRampPalette(rev(brewer.pal(8, color.panel)))
      }else{
        myColorPalette <- colorRampPalette(brewer.pal(8, color.panel))
      }
      gene.median.expression <- aggregate(t(expression.matrix[genes,]), list(Group = groupid), mean)
      gene.median.expression[, 2:(length(genes) + 1)] <- apply(gene.median.expression[, 2:(length(genes) + 1)], 2, function(x){
        (x - mean(x))/sd(x)
      })
      heatmap.data <- reshape2::melt(gene.median.expression, id = "Group")
      heatmap.data$Group <- as.factor(heatmap.data$Group)
      heatmap.data$Group <- factor(heatmap.data$Group, levels = rev(levels(heatmap.data$Group)))
      max_value <- round(max(abs(heatmap.data$value)) + 0.05, 1)
      scale_breaks <- round(max_value * c(-4/5, -2/5, 0, 2/5, 4/5), 1)
      p_heatmap <- ggplot(heatmap.data, aes(x = factor(variable), y = Group)) +
        geom_tile(aes(fill = value), colour = "white") + 
        theme_bw() + 
        theme(panel.border = element_blank(),
              legend.title = element_blank(),
              axis.line = element_blank(),
              #axis.text.y = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                         hjust = 1, size = font.size * 0.8),
              axis.title = element_blank(),
              axis.ticks = element_blank()) +
        scale_fill_gradientn(limits = c(-max_value, max_value), breaks = scale_breaks, expand = c(0,0), colours = myColorPalette(100), guide = guide_colourbar(label.theme = element_text(size = font.size * 0.5))) +
        scale_x_discrete(breaks = unique(heatmap.data$variable), labels = unique(heatmap.data$variable)) + 
        scale_y_discrete(expand = c(0,0))
      if (coord.flip) {
        p_heatmap <- p_heatmap + coord_flip()
      }
      if (!only.heatmap) {p_heatmap <- p_heatmap + theme(axis.text.y = element_blank())}
      gt = ggplotGrob(p_heatmap)
      leg = gtable::gtable_filter(gt, "guide-box")
      leg[[1]][[1]][[1]][[1]][[1]][[2]]$height = unit(1, "npc")
      pos = unit.c(unit(0.01,"npc"), unit(.25, "npc"), unit(.5, "npc"), unit(.75, "npc"), unit(.99, "npc"))
      leg[[1]][[1]][[1]][[1]][[1]][[5]]$y0 = pos
      leg[[1]][[1]][[1]][[1]][[1]][[5]]$y1 = pos
      leg[[1]][[1]][[1]][[1]]$heights = unit.c(rep(unit(0, "mm"), 3),
                                               unit(1, "npc"),
                                               unit(0, "mm"))
      leg[[1]][[1]]$heights[[3]] = sum(rep(unit(0, "mm"), 3),
                                       unit(1, "npc"),
                                       unit(0, "mm"))
      if(strsplit(as.character(packageVersion("ggplot2")),"[.]")[[1]][1] == "2"){
        # ggplot 2.0
        leg[[1]][[1]][[1]][[1]][[1]][[3]]$y = pos
        p_heatmap = gtable::gtable_add_grob(gt, leg, t = 6, l = 8)
      }else if(strsplit(as.character(packageVersion("ggplot2")),"[.]")[[1]][1] == "3"){
        # ggplot 3.0
        leg[[1]][[1]][[1]][[1]][[1]][[3]]$children[[1]]$y = pos
        p_heatmap = gtable::gtable_add_grob(gt, leg, t = 7, l = 9)
      }
      # Final plot
      if (only.heatmap) {
        grid.newpage()
        return(grid.draw(p_heatmap))
      }
      g1 <- ggplotGrob(p_boxplot)
      g2 <- p_heatmap
      fg1 <- egg::gtable_frame(g1, height = unit(10, "null"))
      fg2 <- egg::gtable_frame(g2, height = unit(40, "null"))
      p <- egg::gtable_frame(gridExtra::gtable_cbind(fg1, fg2))
      grid.newpage()
      grid.draw(p)
    } else{
      ggplot(data.frame()) + 
        ggtitle("There should be more than two groups") + 
        theme_bw() + 
        theme(plot.title = element_text(
          size = 26,
          face = "bold",
          hjust = 0.5)
        )
    }
  }else{
    ggplot(data.frame()) +
      ggtitle("There should be more than two genes") +
      theme_bw() +
      theme(plot.title = element_text(
        size = 26,
        face = "bold",
        hjust = 0.5
      ))
  }
}
