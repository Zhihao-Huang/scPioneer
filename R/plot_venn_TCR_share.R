#' Transparent colors
#' Inspired by Mark Gardener 2015, www.dataanalytics.org.uk
#' 
#' @export
t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  colort <- sapply(color, function(i) {
    ## Get RGB values for named color
    rgb.val <- col2rgb(i)
    ## Make new color using input color as base and alpha set by transparency
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                 max = 255,
                 alpha = (100 - percent) * 255 / 100,
                 names = name)
    ## Save the color
    invisible(t.col)
  })
  names(colort)
}
#' Venn plot of shared TCR
#' 
#' @export
vennplot <- function(meta, group, clonotype, alpha = 0.5,margin = 0.2,
                     cat.pos = NULL,
                     main = NULL, filename = NULL,...) {
  allTCR_list <- lapply(levels(meta[,group]), function(x) {
    unique(meta[,clonotype][meta[,group] == x])
  })
  names(allTCR_list) <- levels(meta[,group])
  num_g <- length(allTCR_list)
  if (is.null(cat.pos)) {
    
  }
  #plot
  p <- VennDiagram::venn.diagram(allTCR_list, filename = filename, main = main,
                    fill = scPalette1(num_g), alpha = alpha, 
                    cat.col = rep('black', num_g), cat.default.pos = 'outer',
                    #cat.pos = cat.pos,
                    col = 'black', cex = 1, fontfamily = 'serif', cat.cex = 1,
                    cat.fontfamily = 'serif', margin = margin,...)
  if (is.null(filename)) {
    return(p)
  }
}
#' Venn plot + bar plot of shared TCR
#' 
#' @export
plot_stat <- function(meta, group, clonotype, alpha = 0.5,margin = 0.2,
                      main = NULL, vennfilename = NULL,...) {
  if (is.null(levels(meta[,group]))) {
    meta[,group] <- factor(meta[,group], levels = unique(meta[,group]))
  }
  venn <- vennplot(meta, group, clonotype, margin = margin,
                   alpha = alpha, main = main,...)
  metau <- unique(meta[,c(group, clonotype)])
  plotdata <- metau %>% group_by_(group, .drop = FALSE) %>% dplyr::summarise(count = n())
  # plotdata[,group] <- factor(plotdata[,group], levels = rev(unique(plotdata[,group])))
  #color <- t_col(scPalette1(nrow(plotdata)), percent = alpha - 0.05)
  gg <- ggplot(plotdata, aes_string(group, 'count', fill = group)) +
    geom_bar(position = 'dodge',stat = 'identity') +
    geom_text(aes(label=count), vjust=0) +
    scale_fill_manual(values = alpha(scPalette1(nrow(plotdata)), alpha)) +
    theme_classic() + guides(fill = F) +
    xlab('') + ylab('Total counts') + ylim(0, max(plotdata$count) * 1.2)
  plist <- list('venn' = venn, 'bar' = gg)
  return(plist)
  #grid.newpage()
  #grid.draw(plist[['venn']])
}