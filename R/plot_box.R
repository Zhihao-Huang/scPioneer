#' boxplot that enables row and col facet.
#' @export
boxplot2 <- function(data = NULL, paths = NULL, colGroup =  NULL, rowGroup = NULL,
                     meltdata = NULL, facet.toward = 'row',
                     do.log = F,strip.text.size = 8,color.use = NULL, jitter = F, cell.order = NULL,
                     add.ave.point = F, add.line = F, line.size = 0.5, add.point = T,
                     mid.color.point = 'auto', legend.position = 'top',legend.name = '',
                     jitter.size = 0.1,  text.x.size = 8, text.y.size = 8,
                     hjust = 0, vjust = 0, legend.title.size = 8, axis.text.x.angle = NULL,
                     do.test = F, test.method = c('wilcox.test','t.test','anova','kruskal.test')[1]) {
  if (is.null(meltdata)) {
    sig_path <- as.data.frame(data[,paths, drop=FALSE])
    if (do.log) {
      sig_path <- log2(sig_path + 1)
    }
    sig_path$colGroup <- data[,colGroup]
    sig_path$rowGroup <- data[,rowGroup]
    matm <- melt(sig_path, id.vars = c('colGroup','rowGroup'), measure.vars = paths,
                 variable.name = "pathway", value.name = "TPM")
  }else{
    matm <- meltdata
  }
  if (!is.null(cell.order)) {
    matm$groups <- factor(matm$colGroup, levels = cell.order)
  }
  color_palette = c("#00A087FF", "#4DBBD5FF", "#E64B35FF", "#3C5488FF",
                    "#F38400", "#A1CAF1", "#BE0032", "#C2B280", "#848482", "#008856",
                    "#E68FAC", "#0067A5", "#604E97", "#F6A600", "#B3446C", "#DCD300",
                    "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26","#848482",
                    "#008856", "#E68FAC", "#0067A5", "#604E97", "#F6A600", "#B3446C",
                    "#DCD300", "#882D17")
  p <- ggplot(matm, aes(x = colGroup, y = TPM, color = colGroup))
  p <- p + geom_boxplot(fill= NA) +  theme_bw()
  if (is.null(color.use) & length(unique(matm$groups)) < 31) {
    p <- p + scale_color_manual(values = color_palette)
  }else{
    p <- p + scale_color_manual(values = color.use)
  }
  legend.position <- 'none'
  if (jitter){
    p <- p + geom_jitter(size = jitter.size)
  }
  if (add.point) {
    p <- p + geom_point(data = matm, aes(x = colGroup, y = TPM, color = colGroup),
                        shape = 15, size = jitter.size)
  }
  if (facet.toward == 'row') {
    axis.text.x.angle <- ifelse(is.null(axis.text.x.angle), 90, axis.text.x.angle)
    p <- p + facet_grid(rows = vars(pathway), scales = "free") +
      theme(legend.position = legend.position, panel.grid.major = element_line(colour = NA),
            panel.grid.minor = element_line(colour = NA)) +
      theme(axis.text.x = element_text(angle = axis.text.x.angle, size = text.x.size, hjust = hjust, vjust = vjust),
            axis.text.y = element_text(size = text.y.size))
  }else if (facet.toward == 'col') {
    axis.text.x.angle <- ifelse(is.null(axis.text.x.angle), 270, axis.text.x.angle)
    p <- p + coord_flip()
    p <- p + facet_grid(cols = vars(pathway), scales = "free") +
      theme(legend.position = legend.position, panel.grid.major = element_line(colour = NA),
            panel.grid.minor = element_line(colour = NA)) +
      theme(axis.text.x = element_text(angle = axis.text.x.angle, size = text.x.size, hjust = hjust, vjust = vjust),
            axis.text.y = element_text(size = text.y.size))
  }else if (facet.toward == 'all') {
    axis.text.x.angle <- ifelse(is.null(axis.text.x.angle), 270, axis.text.x.angle)
    p <- p + facet_grid(rows =vars(pathway), cols = vars(rowGroup), scales = "free") +
      theme(legend.position = legend.position, panel.grid.major = element_line(colour = NA),
            panel.grid.minor = element_line(colour = NA)) +
      theme(axis.text.x = element_text(angle = axis.text.x.angle, size = text.x.size, hjust = hjust, vjust = vjust),
            axis.text.y = element_text(size = text.y.size))
  }
  p <- p + xlab("") + ylab("")
  p <- p + theme(legend.title = element_text(size = legend.title.size),
                 strip.text = element_text(size = strip.text.size),strip.background = element_blank(),
                 strip.placement = 'outside')
  p <- p + guides(color = FALSE)
  if (do.test) {
    group <- unique(matm[,c('colGroup','rowGroup')])
    grouplist <- lapply(unique(group$rowGroup), function(x) group$colGroup[group$rowGroup == x])
    p <- p + ggpubr::stat_compare_means(mapping = aes_string("colGroup"), method = test.method,
                                        comparisons =unique(grouplist),
                                        bracket.size = 0.3,label.y.npc = 0.7,vjust=0.5,
                                        symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.1),
                                                         symbols = c("****", "***", "**", "*", "ns")))
  }
  return(p)
}

#' boxplot
#' @export
boxplot <- function(data = NULL, paths = NULL, ident = NULL, meltdata =  NULL, facet.toward = 'row',  do.log = F,strip.text.size = 8,
                    color.use = NULL, jitter = F, cell.order = NULL,
                    add.ave.point = F, add.line = F, line.size = 0.5, add.point = T,
                    mid.color.point = 'auto', legend.position = 'top',legend.name = '',
                    jitter.size = 0.1,  text.x.size = 8, text.y.size = 8,
                    hjust = 0, vjust = 0, legend.title.size = 8, axis.text.x.angle = NULL,
                    do.test = F,comparisons = NULL,
                    test.method = c('wilcox.test','t.test','anova','kruskal.test')[1]) {
  if (is.null(meltdata)) {
    sig_path <- as.data.frame(data[,paths, drop=FALSE])
    if (do.log) {
      sig_path <- log2(sig_path + 1)
    }
    sig_path$groups <- ident
    matm <- melt(sig_path, id.vars = c("groups"), measure.vars = paths,
                 variable.name = "pathway", value.name = "TPM")
  }else{
    matm <- meltdata
  }
  # attach(matm)
  colnames(matm) <- c("groups", "pathway", "TPM")
  ave <- tapply(matm$TPM, list(matm$groups, matm$pathway), mean)
  ave <- as.data.frame(ave)
  ave$groups <- rownames(ave)
  ave <- melt(ave, id.vars = c("groups"), measure.vars = paths,
              variable.name = "pathway", value.name = "ave",factorsAsStrings = T)
  matm <- suppressMessages(inner_join(matm,ave))
  if (!is.null(cell.order)) {
    matm$groups <- factor(matm$groups, levels = cell.order)
  }
  color_palette = c("#00A087FF", "#4DBBD5FF", "#E64B35FF", "#3C5488FF",
                    "#F38400", "#A1CAF1", "#BE0032", "#C2B280", "#848482", "#008856",
                    "#E68FAC", "#0067A5", "#604E97", "#F6A600", "#B3446C", "#DCD300",
                    "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26","#848482",
                    "#008856", "#E68FAC", "#0067A5", "#604E97", "#F6A600", "#B3446C",
                    "#DCD300", "#882D17")
  p <- ggplot(matm, aes(x = groups, y = TPM, color = groups))
  p <- p + geom_boxplot(fill= NA) +  theme_bw()
  if (is.null(color.use) & length(unique(matm$groups)) < 31) {
    p <- p + scale_color_manual(values = color_palette)
  }else{
    p <- p + scale_color_manual(values = color.use)
  }
  legend.position <- 'none'
  if (jitter){
    p <- p + geom_jitter(size = jitter.size)
  }
  if (add.point) {
    p <- p + geom_point(data = matm, aes(x = groups, y = TPM, color = groups),
                        shape = 15, size = jitter.size)
  }
  if (add.ave.point) {
    p <- p + geom_point(data = matm, aes(x = groups, y = ave, group = 1), color = 'Black',
                        shape = 15, size = 1.5)
  }
  if (add.line) {
    p <- p + geom_line(data = matm, aes(x = groups, y = ave, group = 1), color = 'Black',
                       size = line.size, linetype = "dashed")
  }
  if (facet.toward == 'row') {
    axis.text.x.angle <- ifelse(is.null(axis.text.x.angle), 90, axis.text.x.angle)
    p <- p + facet_grid(rows = vars(pathway), scales = "free") +
      theme(legend.position = legend.position, panel.grid.major = element_line(colour = NA),
            panel.grid.minor = element_line(colour = NA)) +
      theme(axis.text.x = element_text(angle = axis.text.x.angle, size = text.x.size, hjust = hjust, vjust = vjust),
            axis.text.y = element_text(size = text.y.size))
  }else if (facet.toward == 'col') {
    axis.text.x.angle <- ifelse(is.null(axis.text.x.angle), 270, axis.text.x.angle)
    p <- p + coord_flip()
    p <- p + facet_grid(cols = vars(pathway), scales = "free") +
      theme(legend.position = legend.position, panel.grid.major = element_line(colour = NA),
            panel.grid.minor = element_line(colour = NA)) +
      theme(axis.text.x = element_text(angle = axis.text.x.angle, size = text.x.size, hjust = hjust, vjust = vjust),
            axis.text.y = element_text(size = text.y.size))
  }
  p <- p + xlab("") + ylab("")
  p <- p + theme(legend.title = element_text(size = legend.title.size),
                 strip.text = element_text(size = strip.text.size),strip.background = element_blank(),
                 strip.placement = 'outside')
  p <- p + guides(color = FALSE)
  if (do.test) {
    if(!is.null(comparisons)){
      p <- p + ggpubr::stat_compare_means(mapping = aes_string("State"), method = test.method,
                                          comparisons = comparisons,
                                          bracket.size = 0.3,label.y.npc = 0.7,vjust=0.5,
                                          symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.1),
                                                           symbols = c("****", "***", "**", "*", "ns")))
    }else{
      p <- p + ggpubr::stat_compare_means(mapping = aes_string("State"), method = test.method,
                                          comparisons = list(c(unique(ident)[1],unique(ident)[2])),
                                          bracket.size = 0.3,label.y.npc = 0.7,vjust=0.5,
                                          symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.1),
                                                           symbols = c("****", "***", "**", "*", "ns")))
    }
    
  }
  return(p)
}

#' Boxplot of fraction of cell number.
#' 
#' @examples 
#' metaf <- meta200[,c('Sample','Group','Annotation')]
#' plotbox(metaf)
#' 
#' @export
plotbox <- function (metaf, color = NULL,
                     axis.line.size = 0.1,
                     box.line.size = 0.1,
                     outlier.size = 0.1) {
  colnames(metaf) <- c("Sample", "Group", "Celltype")
  counts <- table(metaf[, c("Sample", "Celltype")])
  frac <- apply(counts, 1, function(x) x/sum(x))
  gs <- metaf[, c("Sample", "Group")]
  gs <- gs[!duplicated(gs$Sample), ]
  rownames(gs) <- gs$Sample
  gs <- gs[colnames(frac), ]
  fract <- as.data.frame(t(frac))
  fract$Group <- gs$Group
  fract <- reshape2::melt(fract, id.vars = c("Group"), measure.vars = rownames(frac))
  colnames(fract) <- c("Group", "Celltype", "Frac")
  fracb <- as.data.frame(t(frac))
  fracb$Group <- gs$Group
  fracb$Sample <- gs$Sample
  plotdata <- melt(fracb, id.vars = c("Sample", "Group"), measure.vars = rownames(frac))
  colnames(plotdata) <- c("Sample", "Group", "Celltype", "Frac")
  if (is.null(color)) {
    color <- scPioneer::scPalette1(length(unique(plotdata$Group)))
  }
  p <- ggplot(plotdata, aes(Celltype, Frac, fill = Group)) + 
    geom_boxplot(outlier.size = outlier.size, lwd = box.line.size) + 
    theme_classic() + 
    xlab("Cell types") + 
    ylab("Fraction of cell number") + 
    theme(axis.text.x = element_text(size = 8, color = "black", 
                                     vjust = 1, hjust = 1, angle = 40),
          panel.grid.major = element_line(colour = NA), 
          panel.grid.minor = element_line(colour = NA),
          panel.border = element_blank(),
          line = element_line(size = axis.line.size)) + 
    scale_fill_manual(values = color)
  return(p)
}

#' Boxplot of fraction of cell number.
#' 
#' @examples 
#' metaf <- meta200[,c('Sample','Group','Annotation')]
#' plotbox_facet(metaf)
#' 
#' @export
plotbox_facet <- function (metaf, color = NULL,
                           test = c("t","wilcox"),
                           alternative = c("two.sided", "great", "less"),
                           adjust.method = c("BH","holm", "hochberg", "hommel", 
                                             "bonferroni","BY", "fdr", "none"),
                           labelp = c('p.adj', 'p.adj.signif'),
                           scale.frac = F,
                           ref.group = NULL,
                           ylim = NULL,
                           label.expansion = c(0.05, 0.1),
                           y.position = NULL,
                           order.group = NULL,
                           return.data = F,
                           ...
                           ) {
  test <- match.arg(arg = NULL, choices = test)
  alternative <- match.arg(arg = NULL, choices = alternative)
  adjust.method <- match.arg(arg = NULL, choices = adjust.method)
  labelp <- match.arg(arg = NULL, choices = labelp)
  colnames(metaf) <- c("Sample", "Group", "Celltype")
  counts <- table(metaf[, c("Sample", "Celltype")])
  frac <- apply(counts, 1, function(x) x/sum(x))
  if (scale.frac) {
    frac <- t(apply(frac,1,scale))
    colnames(frac) <- rownames(counts)
    ylabel <- 'Scaled fraction'
  }
  gs <- metaf[, c("Sample", "Group")]
  gs <- gs[!duplicated(gs$Sample), ]
  rownames(gs) <- gs$Sample
  gs <- gs[colnames(frac), ]
  fract <- as.data.frame(t(frac))
  fract$Group <- gs$Group
  fract <- reshape2::melt(fract, id.vars = c("Group"), measure.vars = rownames(frac))
  colnames(fract) <- c("Group", "Celltype", "Frac")
  fracb <- as.data.frame(t(frac))
  fracb$Group <- gs$Group
  fracb$Sample <- gs$Sample
  plotdata <- melt(fracb, id.vars = c("Sample", "Group"), measure.vars = rownames(frac))
  colnames(plotdata) <- c("Sample", "Group", "Celltype", "Frac")
  if (is.null(color)) {
    color <- scPioneer::scPalette1(length(unique(plotdata$Group)))
  }
  plotdata$Group <- as.vector(plotdata$Group)
  if (test == 't') {
    stat.test <- plotdata %>%
      group_by(Celltype) %>%
      rstatix::t_test(Frac ~ Group, alternative = alternative, 
                      ref.group = ref.group) %>%
      rstatix::adjust_pvalue(method = adjust.method) %>%
      rstatix::add_significance('p.adj') %>% rstatix::add_y_position()
  }else{
    stat.test <- plotdata %>%
      group_by(Celltype) %>%
      rstatix::wilcox_test(Frac ~ Group, alternative = alternative,
                           ref.group = ref.group) %>%
      rstatix::adjust_pvalue(method = adjust.method) %>%
      rstatix::add_significance('p.adj') %>% rstatix::add_y_position()
  }
  if (!is.null(order.group)) {
    plotdata$Group <- factor(plotdata$Group, levels = order.group)
  }
  if (is.null(ylim)) {
    maxnum <- min(max(plotdata$Frac) + 0.1, 1)
    ylim <- c(0,maxnum)
  }
  
  # Visualization
  p <- ggpubr::ggboxplot(plotdata, x = "Group", y = "Frac",
    color = "Group", palette = "jco", facet.by = "Celltype",
    ylim = ylim, ...) + ylab('Fraction') + scale_color_manual(values = color)
  if (!is.null(y.position)) {
    p <- p + ggpubr::stat_pvalue_manual(stat.test, label = labelp,
                                        y.position = y.position)
  }else{
    p <- p + ggpubr::stat_pvalue_manual(stat.test, label = "p.adj.signif",
                                        tip.length = 0.01) + 
      scale_y_continuous(expand = expansion(mult = label.expansion))
  }
  if (return.data) {
    datalist <- list(plotdata=plotdata, stat.test=stat.test, p = p)
    return(datalist)
  }else{
    return(p)
  }
}


#' Boxplot of fraction of cell number vs total number of cells based on group.
#' 
#' @examples 
#' metaf <- meta200[,c('Sample','Group','Annotation','maingroups), 
#' total.number.based = 'maingroups']
#' plotbox(metaf)
#' 
#' @export
plotbox2 <- function (metaf, 
                     total.number.based = 'Sample',
                     color = NULL,
                     axis.line.size = 0.1,
                     box.line.size = 0.1,
                     outlier.size = 0.1) {
  #check colnames
  colname <- unique(c("Sample", "Group", "Celltype", total.number.based))
  pos <- colname %in% colnames(metaf)
  if (any(!pos)) stop(paste0('Colname not in data: ', paste(colname[pos], collapse = ',')))
  metaf <- metaf[, colname]
  if (!is.factor(metaf$Group)) metaf$Group <- factor(metaf$Group, 
                                                     levels = sort(unique(metaf$Group)))
  # calculate fraction
  if (total.number.based == 'Sample') {
    frac <- metaf %>% group_by_('Group','Sample', 'Celltype') %>% 
      summarise(total_num = n(), .groups = 'drop') %>% 
      group_by_('Group','Sample') %>% 
      mutate(base_num = sum(total_num)) %>% mutate(Frac = total_num/base_num)
  }else if (total.number.based == 'Group') {
      frac <- metaf %>% group_by_('Group', 'Sample', 'Celltype') %>% 
        summarise(total_num = n(), .groups = 'drop') %>% 
        group_by_('Group') %>% 
        mutate(base_num = sum(total_num)) %>% mutate(Frac = total_num/base_num)
  }else{
    frac <- metaf %>% group_by_('Group','Sample', total.number.based, 'Celltype') %>% 
      summarise(total_num = n(), .groups = 'drop') %>% 
      group_by_('Group','Sample', total.number.based) %>% 
      mutate(base_num = sum(total_num)) %>% mutate(Frac = total_num/base_num)
  }
  # add group to fraction df.
  #if (!'Group' %in% colnames(frac)) {
  #  gs <- unique(metaf[, c("Sample", "Group")])
  #  frac <- left_join(frac, gs, by = 'Sample')
  #}
  frac <- frac[,c("Sample", "Group", "Celltype", "Frac")]
  # create data frame of all group and celltypes to keep the possible 0 value.
  Celltypedf <- data.frame(Celltype = unique(metaf$Celltype))
  gs <- unique(metaf[, c("Sample", "Group")])
  plotdata <- lapply(1:nrow(gs), function(x) {
    df <- Celltypedf
    df$Sample <- gs$Sample[x]
    df$Group <- gs$Group[x]
    df
  })
  plotdata <- do.call(rbind, plotdata)
  plotdata <- left_join(plotdata, frac, by = c("Sample", "Group", "Celltype"))
  plotdata$Frac[is.na(plotdata$Frac)] <- 0
  # convert the group to factor format to keep the possible 0 value.
  #plotdata$Group <- factor(plotdata$Group, levels = levels(metaf$Group))
  if (is.null(color)) {
    color <- scPioneer::scPalette1(length(unique(plotdata$Group)))
  }
  p <- ggplot(plotdata, aes(Celltype, Frac, fill = Group)) + 
    geom_boxplot(outlier.size = outlier.size, lwd = box.line.size) + 
    theme_classic() + 
    xlab("Cell types") + 
    ylab("Fraction of cell number") + 
    theme(axis.text.x = element_text(size = 8, color = "black", 
                                     vjust = 1, hjust = 1, angle = 40),
          panel.grid.major = element_line(colour = NA), 
          panel.grid.minor = element_line(colour = NA),
          panel.border = element_blank(),
          line = element_line(size = axis.line.size)) + 
    scale_fill_manual(values = color)
  return(p)
}

#' Boxplot of fraction of cell number. Scale fraction in each cell type.
#' 
#' @examples 
#' metaf <- meta200[,c('Sample','Group','Annotation')]
#' plotbo3x(metaf)
#' 
#' @export
plotbox3 <- function (metaf, color = NULL, scale.frac = F,
                      ylabel = "Fraction of cell number",
                     axis.line.size = 0.1,
                     box.line.size = 0.1,
                     outlier.size = 0.1) {
  colnames(metaf) <- c("Sample", "Group", "Celltype")
  counts <- table(metaf[, c("Sample", "Celltype")])
  frac <- apply(counts, 1, function(x) x/sum(x))
  if (scale.frac) {
    frac <- t(apply(frac,1,scale))
    colnames(frac) <- rownames(counts)
    ylabel <- 'Scaled fraction'
  }
  gs <- metaf[, c("Sample", "Group")]
  gs <- gs[!duplicated(gs$Sample), ]
  rownames(gs) <- gs$Sample
  gs <- gs[colnames(frac), ]
  fract <- as.data.frame(t(frac))
  fract$Group <- gs$Group
  fract <- reshape2::melt(fract, id.vars = c("Group"), measure.vars = rownames(frac))
  colnames(fract) <- c("Group", "Celltype", "Frac")
  fracb <- as.data.frame(t(frac))
  fracb$Group <- gs$Group
  fracb$Sample <- gs$Sample
  plotdata <- melt(fracb, id.vars = c("Sample", "Group"), measure.vars = rownames(frac))
  colnames(plotdata) <- c("Sample", "Group", "Celltype", "Frac")
  if (is.null(color)) {
    color <- scPioneer::scPalette1(length(unique(plotdata$Group)))
  }
  p <- ggplot(plotdata, aes(Celltype, Frac, fill = Group)) + 
    geom_boxplot(outlier.size = outlier.size, lwd = box.line.size) + 
    theme_classic() + 
    xlab("Cell types") + 
    ylab(ylabel) + 
    theme(axis.text.x = element_text(size = 8, color = "black", 
                                     vjust = 1, hjust = 1, angle = 40),
          panel.grid.major = element_line(colour = NA), 
          panel.grid.minor = element_line(colour = NA),
          panel.border = element_blank(),
          line = element_line(size = axis.line.size)) + 
    scale_fill_manual(values = color)
  return(p)
}
