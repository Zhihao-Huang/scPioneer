#' Function for paired pie plot
#' 
#' @param meta Data frame including group and element scores.
#' @param cellID Select rows by rownames
#' @param sctypes colnames that are used for pie elements.
#' @param group.colname colname of group.
#' @param group.names 2 group names included in group colname.
#' @param min.pct minmum percentage that determine to change names of elements to 'Others'.
#' @param test.use Whether to assign asterisks to legend text.
#' @param change.anno.df A data frame to adjust element names.
#' 
#' @export
pie_plot_group <- function(meta, cellID, sctypes, group.colname, 
                           group.names = c('intact gut', 'leaky gut'),
                           min.pct = 0.5,
                           test.use = c("none", "t", "wilcox"), 
                           adjust.method = c("BH","holm", "hochberg",
                                             "hommel", "bonferroni",
                                            "BY", "fdr",
                                             "none"), 
                           plabel = c("p.adj.signif", "p.adj", "p"),
                           colors = NULL, legend.title = NULL,
                           change.anno.df = NULL, cell.order = NULL) {
  test.use <- match.arg(arg = test.use, choices = c("none", "t", "wilcox"))
  adjust.method <- match.arg(arg = adjust.method, 
                             choices = c("BH","holm", "hochberg",
                                         "hommel", "bonferroni",
                                          "BY", "fdr",
                                         "none"))
  plabel <- match.arg(arg = plabel, choices = c("p.adj.signif", "p.adj", "p"))
  df <- meta[cellID,c('Group',sctypes)]
  dfm <- melt(df, id.vars = group.colname)
  colnames(dfm) <- c('Group','Celltypes','Score')
  if (length(group.names) != 2) {stop('Only 2 groups allowed.')}
  dfm <- dfm[dfm$Group %in% group.names, ]
  dfm <- dfm[order(dfm$Score, decreasing = T),]
  if (!is.null(change.anno.df)) {
    dfm <- left_join(dfm, change.anno.df, by = 'Celltypes')
    dfm <- dfm[,c('Group',colnames(change.anno.df)[2],'Score')]
    colnames(dfm) <- c('Group','Celltypes','Score')
  }
  
  dfm$Celltypes <- as.vector(dfm$Celltypes)
  
  if (is.null(colors)) {
    colors <- c('lightgrey',scPalette2(length(unique(dfm$Celltypes))))
    names(colors) <- unique(c('Others',dfm$Celltypes))
  }
  
  plotdata_in <- dfm[dfm$Group == group.names[1],]
  plotdata_in <- plotdata_in %>% group_by(Group, Celltypes) %>% 
    summarise(total_s = sum(Score)) %>%
    mutate(percentage = total_s / sum(total_s) * 100) 
  
  plotdata_in <- plotdata_in[order(plotdata_in$percentage, decreasing = T),]
  
  other_pos <- plotdata_in$percentage < min.pct
  plotdata_in$Celltypes <- as.vector(plotdata_in$Celltypes)
  plotdata_in$percentage <- round(plotdata_in$percentage, 2)
  if (any(other_pos)) {
    plotdata_in$Celltypes[other_pos] <- 'Others'
    alltype_used <- unique(as.vector(plotdata_in$Celltypes[!other_pos]))
    show.levels <- c(alltype_used, 'Others')
  }else{
    alltype_used <- unique(as.vector(plotdata_in$Celltypes))
    show.levels <- alltype_used
  }
  # set 'Others' at the end of levels
  if ('Others' %in% show.levels) show.levels <- rev(unique(c('Others',rev(show.levels))))
  
  plotdata_in$Celltypes <- factor(plotdata_in$Celltypes, levels = show.levels)
  
  
  plotdata_in <- plotdata_in %>% group_by(Group, Celltypes) %>% 
    summarise(pct = sum(percentage)) %>%
    mutate(lab.pos = cumsum(pct)+.5*pct)
  plotdata_in$lab.pos[plotdata_in$pct != max(plotdata_in$pct)] <- NA
  
  #
  
  plotdata_le <- dfm[dfm$Group == group.names[2],]
  plotdata_le <- plotdata_le %>% group_by(Group, Celltypes) %>% 
    #summarise(percentage = mean(Score) * 100) 
    summarise(total_s = sum(Score)) %>%
    mutate(percentage = total_s / sum(total_s) * 100)
  
  plotdata_le$percentage <- round(plotdata_le$percentage, 2)
  other_pos <- !plotdata_le$Celltypes %in% alltype_used
  if (any(other_pos)) {
    plotdata_le$Celltypes[other_pos] <- 'Others'
    alltype_used <- alltype_used[alltype_used %in% plotdata_le$Celltypes]
    show.levels <- unique(c(alltype_used, 'Others'))
  }else{
    show.levels <- alltype_used
  }
  # set 'Others' at the end of levels
  if ('Others' %in% show.levels) show.levels <- rev(unique(c('Others',rev(show.levels))))
  
  plotdata_le$Celltypes <- factor(plotdata_le$Celltypes, levels = show.levels)
  plotdata_le <- plotdata_le[order(plotdata_le$Celltypes),]
  plotdata_le <- plotdata_le %>% group_by(Group, Celltypes) %>% 
    summarise(pct = sum(percentage)) %>%
    mutate(lab.pos = cumsum(pct)+.5*pct)
  plotdata_le$lab.pos[plotdata_le$pct != max(plotdata_le$pct)] <- NA
  
  # do test or not.
  if (test.use != 'none') {
    dfm$Group <- factor(dfm$Group, levels = group.names)
    if (test.use == 't') {
      stat.test <- dfm %>%
        group_by(Celltypes) %>%
        t_test(Score ~ Group) %>%
        adjust_pvalue(method = adjust.method) %>%
        add_significance( "p.adj")
    }
    if (test.use == 'wilcox') {
      stat.test <- dfm %>%
        group_by(Celltypes) %>%
        wilcox_test(Score ~ Group) %>%
        adjust_pvalue(method = adjust.method) %>%
        add_significance( "p.adj")
    }
    statdf <- stat.test[,c('Celltypes', plabel)]
    colnames(statdf) <- c('Celltypes', 'plabel')
    statdf$plabel <- gsub('ns','',statdf$plabel)
    plotdata_le <- left_join(plotdata_le, statdf, 
                             by = 'Celltypes')
    plotdata_le$Celltypes <- as.vector(plotdata_le$Celltypes)
    plotdata_le$Celltypes <- paste0(plotdata_le$Celltypes, '_', plotdata_in$pct,
                                    '% vs ', plotdata_le$pct,'%  ', plotdata_le$plabel) 
  }else{
    plotdata_le$Celltypes <- as.vector(plotdata_le$Celltypes)
    plotdata_le$Celltypes <- paste0(plotdata_le$Celltypes, '_', plotdata_in$pct,
                                    '% vs ', plotdata_le$pct,'%') 
  }

  plotdata_le$Celltypes <- factor(plotdata_le$Celltypes, levels = plotdata_le$Celltypes)
  
  if (!is.null(cell.order)) {
    show.levels <- cell.order
  }
  colordf <- as.data.frame(colors)
  rownames(colordf) <- names(colors)
  color.use <- colordf[show.levels,'colors']
  names(color.use) <- show.levels
  p_in <- ggplot(data = plotdata_in, 
                 aes(x = 2, y = pct, fill = Celltypes))+
    geom_bar(stat = "identity")+
    coord_polar("y", start = 200) +
    #geom_text(aes(y = lab.pos, label = paste(pct,"%", sep = "")), col = "white") +
    theme_void() +
    scale_fill_manual(values = color.use) +
    xlim(.2,2.5) + 
    ggtitle(group.names[1]) + theme(plot.title = element_text(hjust = 0.5),
                                    legend.text = element_text(size=12),
                                    legend.key.size = unit(2, 'cm'))
  
  names(color.use) <- plotdata_le$Celltypes
  p_le <- ggplot(data = plotdata_le, 
                 aes(x = 2, y = pct, fill = Celltypes))+
    geom_bar(stat = "identity")+
    coord_polar("y", start = 200) +
    #geom_text(aes(y = lab.pos, label = paste(pct,"%", sep = "")), col = "white") +
    theme_void() +
    scale_fill_manual(values = color.use)+
    xlim(.2,2.5) + 
    ggtitle(group.names[2]) + theme(plot.title = element_text(hjust = 0.5))
  p_in <- p_in + theme(legend.position = "none")
  if (!is.null(legend.title)) {
    p_le <- p_le + guides(fill=guide_legend(title = legend.title))
  }
  return(p_in + p_le)
}
