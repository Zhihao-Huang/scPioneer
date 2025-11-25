
#' Define M1/M2 phenotypes of Myeloid cells as mean expression of gene signatures, and plot data.
#' M1/M2 markers were from DOI: 10.1016/j.cell.2020.03.048 and DOI: 10.1016/j.cell.2018.05.060 
#' 
#' @param object A Seurat object
#' @param Macrophage.names A vector of all the Macrophage sub-type names.
#' @param return.data Return statistic data.
#' @param plot.type The type of plot. Only 'violin', 'violin merge','box' or 'correlation' is provided.
#' 
#' @return A ggplot object of violin, box or heatmap.
#' 
#' @export
M1M2_score <- function (object, Macrophage.names, return.data = F, 
                        M1M2list = NULL, 
                        M1M2list_Cell2018 = F ,
                        AddModuleScore = F,
                        plot.type = c("violin", 'violin merge', 'box', 'correlation'),
                        celltype.order = NULL,
                        order.by.M1M2 = c('none','M1','M2','diff'),
                        correlation_centre = T, remove.axis.x = F,
                        legend.title = 'Average expression', legend.position = 'right',
                        legend.title.size = 8, 
                        strip.position = 'bottom', nrow = 1,
                        strip.text.size = 8, strip.text.angle = 45, strip.text.hjust = 1, 
                        strip.text.vjust = 1, panel.spacing = 0.5, clip = "on", 
                        border.color = "black", border.size = 0.5) 
{
  plot.type <- match.arg(arg =NULL, choices = plot.type)
  order.by.M1M2 <- match.arg(arg =NULL, choices = order.by.M1M2)
  if (is.null(M1M2list)) {
    Mlist <- list()
    Mlist[["M1"]] <- c("CCL5", "CCR7", "CD40", "CD86", "CXCL9", 
                       "CXCL10", "CXCL11", "IDO1", "IL1A", "IL1B", "IL6", "IRF1", 
                       "IRF5", "KYNU", "KIAA0754", "LILRA3")
    Mlist[["M2"]] <- c("CCL4", "CCL13", "CCL18", "CCL20", "CCL22", 
                       "CD276", "CLEC7A", "CTSA", "CTSB", "CTSC", "CTSD", "FN1", 
                       "IL4R", "IRF4", "LYVE1", "MMP9", "MMP14", "MMP19", "MSR1", 
                       "TGFB1", "TGFB2", "TGFB3", "TNFSF8", "TNFSF12", "VEGFA", 
                       "VEGFB", "VEGFC", "GSTT1")
  }else if (M1M2list_Cell2018) {
    # Refer to doi: 10.1016/j.cell.2018.03.034, table 4
    Mlist <- list()
    Mlist[["M1"]] <- c('IL12B','IL23A','TNF','IL6','CD86',
                       ## 'HLA-DRA','HLA-DRB1','HLA-DRB3','HLA-DRB5','MARCO','PTPRC','CD68','CD14','CLEC10A','CSF1R'
                       ## CXCR10
                       'IL1B','NOS2','FCGR3A','CD80','IL23',
                       'CXCL9','CXCL10','CXCL11','IL1A','IL1B','CCL5','IRF5','IRF1',
                       'CD40','IDO1','KYNU','CCR7'
                       )
    Mlist[["M2"]] <- c('ARG1','ARG2','IL10','FCGR2A','CD163','FCER2','CD200R1',
                       ## 'HLA-DRA','HLA-DRB1','HLA-DRB3','HLA-DRB5','MARCO','PTPRC','CD68','CD14','CLEC10A','CSF1R'
                       'PDCD1LG2','CD274','MRC1','IL1RN','IL1R2','IL4R','CCL4',
                       'CCL13','CCL20','CCL17','CCL18','CCL22','CCL24','LYVE1',
                       'VEGFA','VEGFB','VEGFC','VEGFD','EGF','CTSA','CTSB','CTSC','CTSD',
                       'TGFB1','TGFB2','TGFB3','MMP14','MMP19','MMP9','CLEC7A',
                       'WNT7B','FASL','TNFSF12','TNFSF8','CD276','VTCN1','MSR1',
                       'FN1','IRF4')
  }else{
    Mlist <- M1M2list
  }
  
  MACpos <- Idents(object) %in% Macrophage.names
  plist <- list()
  if (return.data) {
    return(data_score(object, Mlist, MACpos))
  }
  if (plot.type == "violin") {
    for (st in names(Mlist)) {
      genesetname <- paste0("Macrophages ", st, " Score")
      p <- violin_score(object, Mlist[[st]], genesetname,
                        MACpos, AddModuleScore,border.color,border.size)
      plist[[st]] <- p
    }
    return(plist)
  }
  else if (plot.type == 'violin merge') {
    for (st in names(Mlist)) {
      genesetname <- paste0("Macrophages ", st, " Score")
      p <- violin_score(object, Mlist[[st]], genesetname,
                        MACpos, AddModuleScore,border.color, border.size)
      plist[[st]] <- p
    }
    df1 <- plist$M1$data
    df2 <- plist$M2$data
    df1$Type <- 'M1'
    df2$Type <- "M2"
    df <- rbind(df1, df2)
    if (order.by.M1M2 == 'M1') celltype.order <- unique(df1[order(df1$mean_exp,decreasing = T),]$Celltypes)
    if (order.by.M1M2 == 'M2') celltype.order <- unique(df2[order(df2$mean_exp,decreasing = T),]$Celltypes)
    if (order.by.M1M2 == 'diff') {
      unidf1 <- unique(df1[order(df1$Celltypes), c('Celltypes','mean_exp')])
      unidf2 <- unique(df2[order(df2$Celltypes), c('Celltypes','mean_exp')])
      diffexp <- unidf1$mean_exp - unidf2$mean_exp
      celltype.order <- unidf1[order(diffexp, decreasing = T),'Celltypes']
    }
    if (!is.null(celltype.order)) df$Celltypes <- factor(df$Celltypes, levels = celltype.order)
    p <- ggplot(df, aes(x = Type, y = Expression, fill = mean_exp)) + 
      geom_violin(trim = F, scale = "width", color = border.color, linewidth = border.size) + 
      geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.1,
                   color = 'black', linewidth = border.size) +
      theme_classic() + 
      labs(x = "", y = "Macrophages M1/M2 Score") + 
      #facet_grid(cols = vars(Celltypes), scales = "free", switch = strip.text.switch) + 
      facet_wrap(~Celltypes, strip.position = 'bottom',nrow = nrow) +
      theme(legend.position = legend.position, 
            panel.grid.major = element_line(colour = NA),
            panel.grid.minor = element_line(colour = NA)) 
  
  if (remove.axis.x) {
    p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
                   axis.ticks.x = element_blank())
  }
  top_expr <- max(df$Expression)
  median_expr <- (top_expr + min(df$Expression))/3
  mycol <- RColorBrewer::brewer.pal(11, "RdYlBu")
  p <- p + scale_fill_gradient2(name = legend.title, low = mycol[5], mid = mycol[2], 
                                high = mycol[1], midpoint = median_expr)
  p <- p + theme(legend.title = element_text(size = legend.title.size), 
                 strip.text = element_text(size = strip.text.size, angle = strip.text.angle, 
                                           hjust = strip.text.hjust, vjust = strip.text.vjust), 
                 strip.background = element_blank(), strip.placement = "outside", 
                 panel.spacing = unit(panel.spacing, "lines"))
  #p <- p + guides(color = FALSE)
  if (clip == "off") {
    pg <- ggplotGrob(p)
    for (i in which(grepl("strip-", pg$layout$name))) {
      pg$grobs[[i]]$layout$clip <- "off"
    }
    return(pg)
  }
  else {
    return(p)
  }
  }
  else if (plot.type == "box") {
    subm <- subset(object, cells = colnames(object)[MACpos])
    for (st in names(Mlist)) {
      genesetname <- paste0("Macrophages ", st, " Score")
      p <- box_score(object, Mlist[[st]], genesetname, 
                     MACpos, AddModuleScore,  'black', border.size)
      plist[[st]] <- p
    }
    return(plist)
  }
  else if (plot.type == "correlation") {
    return(score_correlation(object, Mlist, MACpos, AddModuleScore,
                             scatter_average = correlation_centre))
  }
}

#' Box plot of gene score.
#' 
#' @examples 
#' cellpos <- pbmc$Annotation %in% c('Naive CD4 T', 'Memory CD4 T','CD8 T','NK')
#' geneset <- c('CCR7','TCF7','SELL')
#' box_score(pbmc, geneset,'Naive lymphocytes',cellpos)
#' geneset <- c('GATA3','PPARG','PTGDR2','IL4','IL5','IL13','CSF1','IL1RL1','IL4R',
#' 'IL17RB','CCR4','CYSLTR1','HPGDS','TESPA1','AHR')
#' box_score(pbmc, geneset,'Type 2 inflammation in Th2 and ILC2',cellpos)
#' 
#' @export
box_score <- function (object, geneset, genesetname, cellpos, AddModuleScore = F,
                       border.color = 'black',border.size = 1) {
  .Deprecated('boxplot_score')
  if (AddModuleScore) {
    # use Seurat::AddModuleScore to calculate score
    object <- Seurat::AddModuleScore(object = object, features = list(geneset),
                                     name = genesetname)
    plot.data <- object@meta.data[cellpos,paste0(gsub(' ','.',genesetname),1), drop = F]
  }else{
    # use average expression as score
    Mgene_in_data <- geneset %in% rownames(object)
    missgene <- geneset[!Mgene_in_data]
    if (any(!Mgene_in_data)) {
      message(paste0(genesetname, " signatures: ", paste(missgene, 
                                                         collapse = ","), " are not in data."))
    }
    mgene <- geneset[Mgene_in_data]
    plot.data <- as.data.frame(t(GetAssayData(object)[mgene,cellpos]))
    plot.data <- as.data.frame(apply(plot.data, 1, mean))
  }
  plot.data$Celltypes <- as.vector(Idents(object))[cellpos]
  colnames(plot.data) <- c("Expression", "Celltypes")
  scorename <- genesetname
  p10 <- ggplot(plot.data, aes(x = Celltypes, y = Expression, 
                               color = Celltypes, fill = Celltypes)) + 
    stat_boxplot(geom = "errorbar", width = 0.5) +
    geom_boxplot(width = 0.5, outlier.size = 0.1,
                 color = border.color, linewidth = border.size)
  p10 <- p10 + theme_classic() + theme(legend.position = "none") + 
    theme(axis.text.x = element_text(size = 10, hjust = 1, 
                                     vjust = 1, angle = 45)) + labs(x = "", y = scorename)
  type_num <- length(unique(plot.data$Celltypes))
  c2 <- rainbow(type_num, alpha = 0.2)
  c3 <- rainbow(type_num, v = 0.7)
  p10 <- p10 + scale_fill_manual(values = c2) + scale_color_manual(values = c3)
  return(p10)
}

#' Violin plot of gene score.
#' 
#' @examples 
#' cellpos <- pbmc$Annotation %in% c('Naive CD4 T', 'Memory CD4 T','CD8 T','NK')
#' geneset <- c('CCR7','TCF7','SELL')
#' violin_score(pbmc, geneset,'Naive lymphocytes',cellpos)
#' 
#' @export
violin_score <- function (object, geneset, genesetname, cellpos, AddModuleScore = F,
                          border.color = 'black', border.size = 1) 
{
  if (AddModuleScore) {
    # use Seurat::AddModuleScore to calculate score
    object <- Seurat::AddModuleScore(object = object, features = list(geneset),
                                     name = genesetname)
    plot.data <- object@meta.data[cellpos,paste0(gsub(' ','.',genesetname),1), drop = F]
  }else{
    # use average expression as score
    Mgene_in_data <- geneset %in% rownames(object)
    missgene <- geneset[!Mgene_in_data]
    if (any(!Mgene_in_data)) {
      message(paste0(genesetname, " signatures: ", paste(missgene, 
                                                         collapse = ","), " are not in data."))
    }
    mgene <- geneset[Mgene_in_data]
    plot.data <- as.data.frame(t(as.matrix(GetAssayData(object)[mgene,cellpos,drop = F])))
    plot.data <- as.data.frame(apply(plot.data, 1, mean))
  }
  plot.data$Celltypes <- as.vector(Idents(object))[cellpos]
  colnames(plot.data) <- c("Expression", "Celltypes")
  ave <- plot.data %>% group_by(Celltypes) %>% summarize(mean_exp = mean(Expression), 
                                                         .groups = "drop")
  plot.data2 <- plyr::join(plot.data, ave, type = "full", by = "Celltypes")
  scorename <- genesetname
  p10 <- ggplot(plot.data2, aes(x = Celltypes, y = Expression, 
                                fill = mean_exp)) + 
    geom_violin(trim = F, scale = "width", color = border.color, linewidth = border.size) + 
    geom_boxplot(width = 0.1, fill = "white", color = border.color, linewidth = border.size)
  p10 <- p10 + theme_classic() + theme(legend.position = "right") + 
    theme(axis.text.x = element_text(size = 10, hjust = 1, 
                                     vjust = 1, angle = 45)) + labs(x = "", y = scorename)
  top_expr <- max(plot.data2$Expression)
  median_expr <- (top_expr + min(plot.data2$Expression))/3
  mycol <- RColorBrewer::brewer.pal(11, "RdYlBu")
  p10 <- p10 + scale_fill_gradient2(low = mycol[5], mid = mycol[2], 
                                    high = mycol[1], midpoint = median_expr)
  return(p10)
}

#' Violin plot of gene score.
#' 
#' @examples 
#' cellpos <- pbmc$Annotation %in% c('Naive CD4 T', 'Memory CD4 T','CD8 T','NK')
#' geneset <- c('CCR7','TCF7','SELL')
#' violin_score(pbmc, geneset,'Naive lymphocytes',cellpos)
#' 
#' @export
violin_score2 <- function (object, geneset, genesetname, 
                           group.by = NULL, 
                           cell.order = NULL,
                           assay = 'RNA',
                           slot = 'data',
                           scale = F,
                           Add.boxplot = T,
                           AddModuleScore = F,
                           test.use = c("none", "t", "wilcox"),
                           adjust.p = T, 
                           plabel = c("p.adj.signif", "p.adj",'p'), 
                           adjust.method = c("holm", "hochberg", "hommel",
                                             "bonferroni", "BH", "BY", 
                                             "fdr",  "none"), 
                           ref.group = NULL,
                           remove.bracket = F,
                           text.pval.size = 3.88,
                           y.position.adjust = 0,
                           color.use = NULL,
                           midpoint = NULL,
                           outlier.size = 0.5,
                           violin.border.color = 'black',
                           box.border.color = 'black',
                           box.border.size = 0.5,
                           return.stat = F,
                           legend.name = NULL,
                           ...) 
{
  test.use <- match.arg(arg = test.use, choices = c("none", "t", "wilcox"))
  adjust.method <- match.arg(arg = adjust.method,
                             choices = c("holm", "hochberg", "hommel",
                                         "bonferroni", "BH", "BY", 
                                         "fdr",  "none"))
  plabel <- match.arg(arg = plabel, 
                      choices = c("p.adj.signif", "p.adj",'p'))
  
  if (AddModuleScore) {
    # use Seurat::AddModuleScore to calculate score
    object <- Seurat::AddModuleScore(object = object, features = list(geneset),
                                     assay = assay,
                                     name = genesetname)
    plot.data <- object@meta.data[,paste0(gsub(' ','.',genesetname),1), drop = F]
  }else{
    # use average expression as score
    Mgene_in_data <- geneset %in% rownames(object)
    missgene <- geneset[!Mgene_in_data]
    if (any(!Mgene_in_data)) {
      message(paste0(genesetname, " signatures: ", paste(missgene, 
                                                         collapse = ","), " are not in data."))
    }
    mgene <- geneset[Mgene_in_data]
    mat <- GetAssayData(object, slot = slot, assay = assay)[mgene,,drop = F]
    mat <- t(as.matrix(mat))
    if (scale) mat <- scale(mat)
    plot.data <- as.data.frame(mat)
    plot.data <- as.data.frame(apply(plot.data, 1, mean))
  }
  if (is.null(group.by)) {
    plot.data$Celltypes <- as.vector(Idents(object))
  }else{
    plot.data$Celltypes <- object@meta.data[,group.by]
  }
  colnames(plot.data) <- c("Expression", "Celltypes")
  
  level <- levels(plot.data$Celltypes)
  if (is.null(level)) level <- unique(plot.data$Celltypes)
  if (!is.null(cell.order)) level <-  cell.order
  
  ave <- plot.data %>% group_by(Celltypes) %>% summarize(mean_exp = mean(Expression), 
                                                         .groups = "drop")
  plot.data2 <- plyr::join(plot.data, ave, type = "full", by = "Celltypes")
  scorename <- genesetname
  plot.data2$Celltypes <- factor(plot.data2$Celltypes, levels = level)
#  p <- ggplot(plot.data2, aes(x = Celltypes, y = Expression, fill = mean_exp)) +
#    geom_violin(trim = F, scale = "width") + 
#    geom_boxplot(width = 0.1, fill = "white")
  p <- ggpubr::ggviolin(data = plot.data2, x = 'Celltypes', y = 'Expression',
                        fill = 'mean_exp', trim = T, scale = "width",
                        color = violin.border.color)
  if (Add.boxplot) {
    p <- p + geom_boxplot(width = 0.1, fill = "white", outlier.size = outlier.size,
                         linewidth = box.border.size, color = box.border.color) 
  }
  p <- p + theme_classic() + theme(legend.position = "right") + 
    theme(axis.text.x = element_text(size = 10, hjust = 1, 
                                     vjust = 1, angle = 45)) + labs(x = "", y = scorename)
  if (is.null(midpoint)) {
    max_avg_exp <- max(plot.data2$mean_exp)
    min_avg_exp <- min(plot.data2$mean_exp)
    midpoint <- (max_avg_exp - min_avg_exp) * 2/3 + min_avg_exp
    #midpoint <- mean(plot.data2$Expression)
  }
  
  if (is.null(color.use)) {
    mycol <- RColorBrewer::brewer.pal(11, "RdYlBu")
  }else{
    mycol <- color.use
  }
  p <- p + scale_fill_gradient2(name = legend.name, low = mycol[5], mid = mycol[2], 
                                    high = mycol[1], midpoint = midpoint)
  # stat
  if (test.use == "t") {
    pairwise.test <- plot.data2 %>% rstatix::t_test(Expression ~ Celltypes, 
                                            p.adjust.method = adjust.method,
                                            ref.group = ref.group)
  }
  if (test.use == "wilcox") {
    pairwise.test <- rstatix::wilcox_test(plot.data2, Expression ~ Celltypes, 
                                          p.adjust.method = adjust.method,
                                          ref.group = ref.group)
  }
  if (test.use != "none") {
    p <- p + ggpubr::stat_pvalue_manual(pairwise.test, label = plabel, 
                                        y.position = seq(1:nrow(pairwise.test)) * max(plot.data2$Expression)/10 + 
                                          max(plot.data2$Expression) + y.position.adjust,
                                        label.size = text.pval.size,
                                        remove.bracket = remove.bracket,...)
  }
  
  if (return.stat) {
    return(pairwise.test)
  }
  return(p)
}
#' Violin plot of gene score. Genes were averaged in each sample.
#' 
#' @export
violin_score3 <- function (object, geneset, genesetname, 
                           group.by = NULL, 
                           sample.group = NULL,
                           cell.order = NULL,
                           scale = F,
                           Add.boxplot = T,
                           Add.jitter = F,
                           jitter.size = 1,
                           AddModuleScore = F,
                           test.use = c("none", "t", "wilcox"),
                           adjust.p = T, 
                           plabel = c("p.adj.signif", "p.adj",'p'), 
                           adjust.method = c("holm", "hochberg", "hommel",
                                             "bonferroni", "BH", "BY", 
                                             "fdr",  "none"), 
                           ref.group = NULL,
                           violin.trim = T,
                           violin.scale = "width",
                           remove.bracket = F,
                           text.pval.size = 3.88,
                           y.position.adjust = 0,
                           color.use = NULL,
                           midpoint = NULL,
                           outlier.size = 0.5,
                           legend.name = NULL,
                           return.stat = F,
                           ...) 
{
  test.use <- match.arg(arg = NULL, choices = test.use)
  adjust.method <- match.arg(arg = NULL, choices = adjust.method)
  plabel <- match.arg(arg = NULL, choices = plabel)
  
  if (AddModuleScore) {
    # use Seurat::AddModuleScore to calculate score
    object <- Seurat::AddModuleScore(object = object, features = list(geneset),
                                     name = genesetname)
    plot.data <- object@meta.data[,paste0(gsub(' ','.',genesetname),1), drop = F]
  }else{
    # use average expression as score
    Mgene_in_data <- geneset %in% rownames(object)
    missgene <- geneset[!Mgene_in_data]
    if (any(!Mgene_in_data)) {
      message(paste0(genesetname, " signatures: ", paste(missgene, 
                                                         collapse = ","), " are not in data."))
    }
    mgene <- geneset[Mgene_in_data]
    plot.data <- as.data.frame(t(as.matrix(object@assays$RNA@data[mgene,,drop = F])))
    
    plot.data <- as.data.frame(apply(plot.data, 1, mean))
    colnames(plot.data) <- 'ave_all_genes'
    if (scale) {
      plot.data <- apply(plot.data, 2, function(x){
        (x - mean(x))/sd(x)
      })
      plot.data <- as.data.frame(plot.data)
      colnames(plot.data) <- 'ave_all_genes'
    }
  }
  if (is.null(group.by)) {
    plot.data$Celltypes <- as.vector(Idents(object))
  }else{
    plot.data$Celltypes <- object@meta.data[,group.by]
  }
  if (!is.null(sample.group)) {
    plot.data$Group <-  object@meta.data[,sample.group]
    plot.data <- plot.data %>% group_by(Celltypes,Group) %>% 
      summarise(Expression = mean(ave_all_genes),.groups = "drop")
    plot.data <- plot.data[,c("Expression", "Celltypes")]
  }
 
  colnames(plot.data) <- c("Expression", "Celltypes")
  
  level <- levels(plot.data$Celltypes)
  if (is.null(level)) level <- unique(plot.data$Celltypes)
  if (!is.null(cell.order)) level <-  cell.order
  
  ave <- plot.data %>% group_by(Celltypes) %>% summarize(mean_exp = mean(Expression), 
                                                         .groups = "drop")
  plot.data2 <- plyr::join(plot.data, ave, type = "full", by = "Celltypes")
  scorename <- genesetname
  
  plot.data2$Celltypes <- factor(plot.data2$Celltypes, levels = level)
  
  #  p <- ggplot(plot.data2, aes(x = Celltypes, y = Expression, fill = mean_exp)) +
  #    geom_violin(trim = F, scale = "width") + 
  #    geom_boxplot(width = 0.1, fill = "white")
  p <- ggpubr::ggviolin(data = plot.data2, x = 'Celltypes', y = 'Expression',
                        fill = 'mean_exp', trim = violin.trim, scale = violin.scale)
  if (Add.boxplot) {
    p <- p + geom_boxplot(width = 0.1, fill = "white", outlier.size = outlier.size) 
  }
  if (Add.jitter) {
    p <- p + geom_jitter(size = jitter.size)
  }
  p <- p + theme_classic() + theme(legend.position = "right") + 
    theme(axis.text.x = element_text(size = 10, hjust = 1, 
                                     vjust = 1, angle = 45)) + labs(x = "", y = scorename)
  if (is.null(midpoint)) {
    max_avg_exp <- max(plot.data2$mean_exp)
    min_avg_exp <- min(plot.data2$mean_exp)
    midpoint <- (max_avg_exp - min_avg_exp) * 2/3 + min_avg_exp
    #midpoint <- mean(plot.data2$Expression)
  }
  
  if (is.null(color.use)) {
    mycol <- RColorBrewer::brewer.pal(11, "RdYlBu")
  }else{
    mycol <- color.use
  }
  p <- p + scale_fill_gradient2(name = legend.name, low = mycol[5], mid = mycol[2], 
                                high = mycol[1], midpoint = midpoint)
  # stat
  if (test.use == "t") {
    pairwise.test <- plot.data2 %>% rstatix::t_test(Expression ~ Celltypes, 
                                                    p.adjust.method = adjust.method,
                                                    ref.group = ref.group)
  }
  if (test.use == "wilcox") {
    pairwise.test <- rstatix::wilcox_test(plot.data2, Expression ~ Celltypes, 
                                          p.adjust.method = adjust.method,
                                          ref.group = ref.group)
  }
  if (test.use != "none") {
    p <- p + ggpubr::stat_pvalue_manual(pairwise.test, label = plabel, 
                                        y.position = seq(1:nrow(pairwise.test)) * max(plot.data2$Expression)/10 + 
                                          max(plot.data2$Expression) + y.position.adjust,
                                        label.size = text.pval.size,
                                        remove.bracket = remove.bracket,...)
  }
  
  if (return.stat) {
    return(pairwise.test)
  }
  return(p)
}

#' Violin plot of gene score. Genes were averaged in each sample.
#' 
#' @export
violin_score4 <- function (object, geneset, genesetname, 
                           group.by = NULL, 
                           facet.by = NULL,
                           Add.boxplot = T,
                           Add.jitter = F,
                           jitter.size = 1,
                           AddModuleScore = F,
                           test.use = c("none", "t", "wilcox"),
                           adjust.p = T, 
                           plabel = c("p.adj.signif", "p.adj",'p'), 
                           adjust.method = c("holm", "hochberg", "hommel",
                                             "bonferroni", "BH", "BY", 
                                             "fdr",  "none"), 
                           ref.group = NULL,
                           violin.trim = T,
                           remove.bracket = F,
                           text.pval.size = 3.88,
                           y.position.adjust = 0,
                           adjust_p_text_gap = 1, 
                           adjust_p_text_height = 1, 
                           ylim = NULL, ylim.adjust = 1.1,
                           color.use = NULL,
                           midpoint = NULL,
                           outlier.size = 0.5,
                           legend.name = NULL,
                           return.stat = F,
                           ...) 
{
  test.use <- match.arg(arg = NULL, choices = test.use)
  adjust.method <- match.arg(arg = NULL, choices = adjust.method)
  plabel <- match.arg(arg = NULL, choices = plabel)
  
  if (AddModuleScore) {
    # use Seurat::AddModuleScore to calculate score
    object <- Seurat::AddModuleScore(object = object, features = list(geneset),
                                     name = genesetname)
    plot.data <- object@meta.data[,paste0(gsub(' ','.',genesetname),1), drop = F]
  }else{
    # use average expression as score
    Mgene_in_data <- geneset %in% rownames(object)
    missgene <- geneset[!Mgene_in_data]
    if (any(!Mgene_in_data)) {
      message(paste0(genesetname, " signatures: ", paste(missgene, 
                                                         collapse = ","), " are not in data."))
    }
    mgene <- geneset[Mgene_in_data]
    plot.data <- as.data.frame(t(as.matrix(object@assays$RNA@data[mgene,,drop = F])))
    plot.data <- as.data.frame(apply(plot.data, 1, mean))
    colnames(plot.data) <- 'ave_all_genes'
  }
  if (is.null(group.by)) {
    plot.data$Celltypes <- as.vector(Idents(object))
  }else{
    plot.data$Celltypes <- object@meta.data[,group.by]
  }
  #return(plot.data)
  scorename <- genesetname
  if (!is.null(facet.by)) {
    plot.data$Group <-  object@meta.data[,facet.by]
    ave <- plot.data %>% group_by(Celltypes,Group) %>% 
      summarise(mean_exp = mean(ave_all_genes),.groups = "drop")
    #ave <- ave[,c("Expression", "Celltypes", "Group")]
    plot.data2 <- plyr::join(plot.data, ave, type = "full", by = c("Celltypes","Group"))
    colnames(plot.data2) <- c("Expression", "Celltypes", "Group",'mean_exp')
    p <- ggpubr::ggviolin(data = plot.data2, x = 'Celltypes', y = 'Expression',facet.by = 'Group',
                          fill = 'mean_exp', trim = violin.trim,...)
    #return(plot.data2)
    # stat
    # remove samples with number < 2, in case of error: "not enough 'y' observations".
    if (test.use != "none") {
      group_num <- plot.data2 %>% group_by(Group, Celltypes) %>% summarize(count = n())
      pos <- group_num$count < 2
      if (any(pos)) {
        message(paste0("Remove pairs ", unique(paste0(group_num$Celltypes[pos],'_',group_num$Group[pos]))," due to not enough 'y' observations. \n"))
        plot.data2 <- plot.data2[!(plot.data2$Celltypes %in% group_num$Celltypes[pos] & plot.data2$Group %in% group_num$Group[pos]), ]
        
      }
      group_num2 <- plot.data2 %>% group_by(Group) %>% summarize(count = length(unique(Celltypes)))
      pos <- group_num2$count < 2
      if (any(pos)) {
        message(paste0("Remove group ", unique(group_num2$Group[pos])," due to not enough 'y' observations. \n"))
        plot.data2 <- plot.data2[!plot.data2$Group %in% group_num2$Group[pos], ]
      }
      }
    if (test.use == "t") {
      pairwise.test <-  plot.data2 %>% group_by(Group) %>% 
        rstatix::t_test(Expression ~ Celltypes,ref.group = ref.group) %>% 
        rstatix::adjust_pvalue(method = adjust.method) %>%
        rstatix::add_significance( "p.adj")
    }
    if (test.use == "wilcox") {
      pairwise.test <-  plot.data2 %>% group_by(Group) %>% 
        rstatix::wilcox_test(Expression ~ Celltypes,ref.group = ref.group) %>% 
        rstatix::adjust_pvalue(method = adjust.method) %>%
        rstatix::add_significance( "p.adj")
    }
    if (test.use != "none") {
      group_num <- pairwise.test %>% group_by(Group) %>% summarize(count = n())
      y.position.adjust <- unlist(as.data.frame(sapply(group_num$count, 
                                                       function(x) 1:x)))
      y.position <- y.position.adjust * max(plot.data2$Expression)/10 * 
        adjust_p_text_gap + max(plot.data2$Expression)
      y.position <- y.position * adjust_p_text_height
      
      if (is.null(ylim)) {
        ylim = c(0, max(y.position) *ylim.adjust)
      }
      p <- p + ggpubr::stat_pvalue_manual(pairwise.test, label = plabel, 
                                          y.position = y.position,
                                          remove.bracket = remove.bracket) + ylim(ylim)
    }
  }else{
    colnames(plot.data) <- c("Expression", "Celltypes")
    ave <- plot.data %>% group_by(Celltypes) %>% summarize(mean_exp = mean(Expression), 
                                                           .groups = "drop")
    plot.data2 <- plyr::join(plot.data, ave, type = "full", by = "Celltypes")
    p <- ggpubr::ggviolin(data = plot.data2, x = 'Celltypes', y = 'Expression',
                          fill = 'mean_exp', trim = violin.trim, scale = "width",...)
    # stat
    if (test.use == "t") {
      pairwise.test <- plot.data2 %>% 
        rstatix::t_test(Expression ~ Celltypes,ref.group = ref.group) %>% 
        rstatix::adjust_pvalue(method = adjust.method) %>%
        rstatix::add_significance( "p.adj")
    }
    if (test.use == "wilcox") {
      pairwise.test <-  plot.data2 %>%
        rstatix::wilcox_test(Expression ~ Celltypes,ref.group = ref.group) %>% 
        rstatix::adjust_pvalue(method = adjust.method) %>%
        rstatix::add_significance( "p.adj")
    }
    if (test.use != "none") {
      p <- p + ggpubr::stat_pvalue_manual(pairwise.test, label = plabel, 
                                          y.position = seq(1:nrow(pairwise.test)) * max(plot.data2$Expression)/10 + 
                                            max(plot.data2$Expression) + y.position.adjust,
                                          label.size = text.pval.size,
                                          remove.bracket = remove.bracket)
    
    }
  }
   
  #  p <- ggplot(plot.data2, aes(x = Celltypes, y = Expression, fill = mean_exp)) +
  #    geom_violin(trim = F, scale = "width") + 
  #    geom_boxplot(width = 0.1, fill = "white")
  if (Add.boxplot) {
    p <- p + geom_boxplot(width = 0.1, fill = "white", outlier.size = outlier.size) 
  }
  if (Add.jitter) {
    p <- p + geom_jitter(size = jitter.size)
  }
  p <- p + theme_classic() + theme(legend.position = "right") + 
    theme(axis.text.x = element_text(size = 10, hjust = 1, 
                                     vjust = 1, angle = 45)) + labs(x = "", y = scorename)
  if (is.null(midpoint)) {
    max_avg_exp <- max(plot.data2$mean_exp)
    min_avg_exp <- min(plot.data2$mean_exp)
    midpoint <- (max_avg_exp - min_avg_exp) * 2/3 + min_avg_exp
    #midpoint <- mean(plot.data2$Expression)
  }
  
  if (is.null(color.use)) {
    mycol <- RColorBrewer::brewer.pal(11, "RdYlBu")
  }else{
    mycol <- color.use
  }
  p <- p + scale_fill_gradient2(name = legend.name, low = mycol[5], mid = mycol[2], 
                                high = mycol[1], midpoint = midpoint)
  
  
  if (return.stat) {
    return(pairwise.test)
  }
  return(p)
}

#' Function for M1M2_score.
#' 
#' Scatter plot to show correlation of 2 geneset.
#' 
#' @export
score_correlation <- function (object, genesetlist, cellpos, AddModuleScore =F,
                               genesetname = NULL,
                               scatter_average = F) 
{
  datalist <- list()
  gname <- names(genesetlist)
  for (st in gname) {
    geneset <- genesetlist[[st]]
    if (is.null(genesetname)) genesetname = 'M1M2'
    if (AddModuleScore) {
      # use Seurat::AddModuleScore to calculate score
      
      object <- Seurat::AddModuleScore(object = object, features = list(geneset),
                                       name = genesetname)
      plot.data <- object@meta.data[cellpos,paste0(gsub(' ','.',genesetname),1), drop = F]
    }else{
      # use average expression as score
      Mgene_in_data <- geneset %in% rownames(object)
      missgene <- geneset[!Mgene_in_data]
      if (any(!Mgene_in_data)) {
        message(paste0(genesetname, " signatures: ", paste(missgene, 
                                                           collapse = ","), " are not in data."))
      }
      mgene <- geneset[Mgene_in_data]
      plot.data <- as.data.frame(t(object@assays$RNA@data[mgene,cellpos]))
      plot.data <- as.data.frame(apply(plot.data, 1, mean))
    }
    plot.data$Celltypes <- as.vector(Idents(object))[cellpos]
    plot.data$Barcodes <- names(Idents(object))[cellpos]
    colnames(plot.data) <- c(st, "Celltypes", "Barcodes")
    colnames(plot.data) <- c("Expression", "Celltypes", "Barcodes")
    ave <- plot.data %>% group_by(Celltypes) %>% summarize(mean_exp = mean(Expression), 
                                                           .groups = "drop")
    plot.data2 <- plyr::join(plot.data, ave, type = "full", 
                             by = "Celltypes")
    colnames(plot.data2)[c(1, 4)] <- paste0(st, "_", colnames(plot.data2)[c(1, 
                                                                            4)])
    datalist[[st]] <- plot.data2
  }
  plotdata <- inner_join(datalist[[gname[1]]], datalist[[gname[2]]], 
                         by = c("Celltypes", "Barcodes"))
  p <- ggplot(plotdata, aes(x = M1_Expression, y = M2_Expression, 
                            color = Celltypes)) + geom_point(size = 0.1)
  p <- p + theme_classic() + theme(legend.position = "right") + 
    guides(colour = guide_legend(override.aes = list(size = 3))) + 
    labs(x = paste0(gname[1], " Signature Expression"), y = paste0(gname[2], 
                                                                   " Signature Expression"))
  p <- p + geom_abline(intercept = 0, slope = 1, color = "blue")
  mycol <- scPalette1(length(unique(plotdata$Celltypes)))
  p <- p + scale_color_manual(values = mycol)
  p <- p + coord_fixed()
  M1sig <- paste0(gname[1], "_Expression")
  M2sig <- paste0(gname[2], "_Expression")
  corr <- round(cor(plotdata[, M1sig], plotdata[, M2sig]), 
                4)
  p <- p + ggtitle(paste0("Correlation = ", corr)) + theme(plot.title = element_text(hjust = 0.5))
  if (scatter_average) {
    ave <- unique(plotdata[, c(2, 4, 6)])
    p <- p + geom_point(data = ave, aes(M1_mean_exp, M2_mean_exp, 
                                        fill = Celltypes), pch = 21, size = 2, color = "black")
    p <- p + scale_fill_manual(values = mycol)
  }
  return(p)
}


#' Box plot of gene score.
#' 
#' @examples 
#' cellt <- c('Naive CD4 T', 'Memory CD4 T','CD8 T','NK')
#' geneset <- c('CCR7','TCF7','SELL')
#' boxplot_score(pbmc, features = geneset,select.celltypes = cellt)
#' # Show number of star of p value.
#' boxplot_score(pbmc, features = geneset,select.celltypes = cellt, 
#' test.use = 'wilicox', adjust.method = 'bonferroni')
#' # Show p value
#' boxplot_score(pbmc, features = geneset,select.celltypes = cellt, 
#' test.use = 'wilicox', adjust.method = 'bonferroni', plabel = 'p.adj')
#' # Return statistic data frame
#' df <- boxplot_score(pbmc, features = geneset,select.celltypes = cellt, 
#' test.use = 'wilicox', adjust.method = 'bonferroni', return.stat = T)
#' 
#' @export
boxplot_score <- function (object, assay = "RNA", features, group.by = NULL,
                           select.celltypes = NULL, 
                           color.use = NULL, 
                           test.use = c("none", "t", "wilicox"), adjust.p = T, 
                           plabel = c("p.adj.signif", "p.adj",'p'), 
                           adjust.method = c("holm", "hochberg", "hommel",
                                             "bonferroni", "BH", "BY", 
                                             "fdr",  "none"), return.stat = F) {
  test.use <- match.arg(arg = NULL, choices = test.use)
  adjust.method <- match.arg(arg = NULL, choices = adjust.method)
  plabel <- match.arg(arg = NULL, choices = plabel)
  genepos <- features %in% c(rownames(object),colnames(object@meta.data))
  if (any(!genepos)) {
    message(paste0("Warning: ", paste(features[!genepos], 
                                      collapse = ", "), " were not in object."))
    features <- features[genepos]
  }
  DefaultAssay(object) <- assay
  data <- Seurat::FetchData(object, features) 
  if (is.null(group.by)) {
    cluster <- Idents(object)
  }
  else {
    cluster <- object@meta.data[, group.by]
  }
  data$celltypes <- cluster
  if (!is.null(select.celltypes)) {
    data$celltypes <- as.vector(data$celltypes)
    pos <- data$celltypes %in% select.celltypes
    if (sum(pos) == 0) {
      stop("All types specified by select.celtypes are not in ident.")
    }
    pos2 <- select.celltypes %in% data$celltypes
    if (sum(pos2) != length(pos2)) {
      message(paste0("Warning: cell types that were not in ident: ", 
                     paste(select.celltypes[!pos2], collapse = ", ")))
    }
    data <- data[pos, ]
    data$celltypes <- factor(data$celltypes, levels = select.celltypes[pos2])
  }
  avemat <- apply(data[-ncol(data)], 1, mean)
  df <- data.frame(data$celltypes, avemat)
  colnames(df) <- c("Celltype", "ave_exp")
  if (class(df$Celltype) == "factor") {
    newlevels <- levels(df$Celltype)[table(df$Celltype) != 
                                       0]
    df$Celltype <- factor(df$Celltype, levels = newlevels)
  }
  if (test.use == "t") {
    pairwise.test <- df %>% rstatix::t_test(ave_exp ~ Celltype, 
                                            p.adjust.method = adjust.method)
  }
  if (test.use == "wilicox") {
    pairwise.test <- rstatix::wilcox_test(df, ave_exp ~ Celltype, 
                                          p.adjust.method = adjust.method)
  }
  p <- ggpubr::ggboxplot(df, x = "Celltype", y = "ave_exp", 
                         color = "Celltype")
  if (test.use != "none") {
    p <- p + ggpubr::stat_pvalue_manual(pairwise.test, label = plabel, 
                                        y.position = seq(1:nrow(pairwise.test)) * max(df$ave_exp)/10 + 
                                          max(df$ave_exp))
  }
  p <- p + ylab("Average expression") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  if (is.null(color.use)) {
    color.use <- scPalette1(length(unique(df$Celltype)))
  }
  p <- p + scale_color_manual(values = color.use)
  if (return.stat) {
    return(pairwise.test)
  }
  p
}

#' Box plot of gene score for each gene.
#' 
#' @examples 
#' cellt <- c('Naive CD4 T', 'Memory CD4 T','CD8 T','NK')
#' subs <- subset(pbmc, Annotation %in% cellt)
#' geneset <- c('CCR7','TCF7','SELL')
#' boxplot_score_facetgene(subs, features = geneset)
#' 
#' @export
boxplot_score_facetgene <- function (object, assay = "RNA", features, 
                                     group.by = NULL, select.group = NULL, 
                                     color.use = NULL, 
                                     test.use = c("none", "t", "wilcox"), 
                                     adjust.method = c("holm", "hochberg",
                                                       "hommel", "bonferroni",
                                                       "BH", "BY", "fdr",
                                                       "none"), 
                                     plabel = c("p.adj.signif", "p.adj", "p"),
                                     ref.group = NULL,
                                     remove.bracket = F,
                                     adjust_p_text_gap = 1, 
                                     adjust_p_text_height = 1, 
                                     ylim = NULL, ylim.adjust = 1.1,
                                     return.stat = F,...) {
  test.use <- match.arg(arg = test.use, choices = c("none", "t", "wilcox"))
  adjust.method <- match.arg(arg = adjust.method, choices = c("holm", "hochberg",
                                                              "hommel", "bonferroni",
                                                              "BH", "BY", "fdr",
                                                              "none"))
  plabel <- match.arg(arg = plabel, choices = c("p.adj.signif", "p.adj", "p"))
  genepos <- features %in% rownames(object)
  if (!is.null(ref.group)) {
    if (!ref.group %in% object@meta.data[, group.by]) ref.group = NULL
  }
  if (any(!genepos)) {
    message(paste0("Warning: ", paste(features[!genepos], 
                                      collapse = ", "), " were not in object."))
    features <- features[genepos]
  }
  data <- as.data.frame(t(as.matrix(GetAssayData(object, assay = assay)[features, 
  ])))
  if (is.null(group.by)) {
    cluster <- Idents(object)
  }
  else {
    cluster <- object@meta.data[, group.by]
  }
  data$celltypes <- cluster
  if (!is.null(select.group)) {
    data$celltypes <- as.vector(data$celltypes)
    pos <- data$celltypes %in% select.group
    if (sum(pos) == 0) {
      stop("All types specified by select.celtypes are not in ident.")
    }
    pos2 <- select.group %in% data$celltypes
    if (sum(pos2) != length(pos2)) {
      message(paste0("Warning: cell types that were not in ident: ", 
                     paste(select.group[!pos2], collapse = ", ")))
    }
    data <- data[pos, ]
    data$celltypes <- factor(data$celltypes, levels = select.group[pos2])
  }
  matm <- reshape2::melt(data)
  colnames(matm) <- c("celltypes", "features", "expression")
  if (test.use == "t") {
    stat.test <- matm %>% group_by(features) %>% 
      rstatix::t_test(expression ~celltypes,
                      p.adjust.method = adjust.method, 
                      ref.group = ref.group) %>%
      add_significance( "p.adj")
  }
  if (test.use == "wilcox") {
    stat.test <- matm %>% group_by(features) %>% 
      rstatix::wilcox_test(expression ~celltypes,
                           p.adjust.method = adjust.method, 
                           ref.group = ref.group) %>%
      add_significance( "p.adj")
  }

  if (test.use != "none") {
    group_num <- stat.test %>% group_by(features) %>% summarize(count = n())
    y.position.adjust <- unlist(as.data.frame(sapply(group_num$count, 
                                                     function(x) 1:x)))
    y.position <- y.position.adjust * max(matm$expression)/10 * 
      adjust_p_text_gap + max(matm$expression)
    y.position <- y.position * adjust_p_text_height
    
    if (is.null(ylim)) {
      ylim = c(0, max(y.position) *ylim.adjust)
    }
    p <- ggpubr::ggboxplot(matm, x = "celltypes", y = "expression",
                           color = "celltypes", 
                           facet.by = "features", ylim = ylim) + xlab("") 
    
    p <- p + ggpubr::stat_pvalue_manual(stat.test, label = plabel, 
                                        y.position = y.position,
                                        remove.bracket = remove.bracket,...)
    
  }else{
    p <- ggpubr::ggboxplot(matm, x = "celltypes", y = "expression",
                           color = "celltypes", 
                           facet.by = "features", ylim = ylim) + xlab("") 
  }
  if (is.null(color.use)) {
    return(p + ggsci::scale_color_jco())
  }
  p <- p + scale_color_manual(values = color.use)
  return(p)
}


#' Box plot of gene score for each group
#' 
#' @examples 
#' cellt <- c('Naive CD4 T', 'Memory CD4 T','CD8 T','NK')
#' pbmc$sample <- c(rep('SampleA',1000), rep('SampleB',1000), rep('SampleC',638))
#' subs <- subset(pbmc, Annotation %in% cellt)
#' geneset <- c('CCR7','TCF7','SELL')
#' boxplot_score_facetgene(subs, features = geneset)
#' 
#' @export
boxplot_score_facetgroup <- function (object, assay = "RNA", features, 
                                     group.by = NULL, select.group = NULL, 
                                     sample.group = NULL,
                                     color.use = NULL, 
                                     test.use = c("none", "t", "wilcox"), 
                                     adjust.method = c("holm", "hochberg",
                                                       "hommel", "bonferroni",
                                                       "BH", "BY", "fdr",
                                                       "none"), 
                                     plabel = c("p.adj.signif", "p.adj", "p"),
                                     ref.group = NULL,
                                     remove.bracket = F,
                                     adjust_p_text_gap = 1, 
                                     adjust_p_text_height = 1, 
                                     ylim = NULL, ylim.adjust = 1.1,
                                     return.stat = F,...) {
  test.use <- match.arg(arg = test.use, choices = c("none", "t", "wilcox"))
  adjust.method <- match.arg(arg = adjust.method, choices = c("holm", "hochberg",
                                                              "hommel", "bonferroni",
                                                              "BH", "BY", "fdr",
                                                              "none"))
  plabel <- match.arg(arg = plabel, choices = c("p.adj.signif", "p.adj", "p"))
  genepos <- features %in% rownames(object)
  if (!is.null(ref.group)) {
    if (!ref.group %in% object@meta.data[, group.by]) ref.group = NULL
  }
  if (any(!genepos)) {
    message(paste0("Warning: ", paste(features[!genepos], 
                                      collapse = ", "), " were not in object."))
    features <- features[genepos]
  }
  data <- as.data.frame(t(as.matrix(GetAssayData(object, assay = assay)[features, 
  ])))
  if (is.null(group.by)) {
    cluster <- Idents(object)
  }
  else {
    cluster <- object@meta.data[, group.by]
  }
  data$celltypes <- cluster
  if (!is.null(select.group)) {
    data$celltypes <- as.vector(data$celltypes)
    pos <- data$celltypes %in% select.group
    if (sum(pos) == 0) {
      stop("All types specified by select.celtypes are not in ident.")
    }
    pos2 <- select.group %in% data$celltypes
    if (sum(pos2) != length(pos2)) {
      message(paste0("Warning: cell types that were not in ident: ", 
                     paste(select.group[!pos2], collapse = ", ")))
    }
    data <- data[pos, ]
    data$celltypes <- factor(data$celltypes, levels = select.group[pos2])
  }
  matm <- reshape2::melt(data)
  colnames(matm) <- c("celltypes", "features", "expression")
  
  if (!is.null(sample.group)) {
    annodf <-  data.frame(celltypes = cluster,
                          sampleg = object@meta.data[,sample.group],
                          stringsAsFactors = F)
    matm <- left_join(matm, annodf)
    matm <- matm %>% group_by(Celltypes,Group) %>% 
      summarise(Expression = mean(expression),.groups = "drop")
    matm <- matm[,c("Expression", "Celltypes")]
  }
  colnames(matm) <- c("Expression", "Celltypes")
  ave <- matm %>% group_by(Celltypes) %>% summarize(mean_exp = mean(Expression), 
                                                    .groups = "drop")
  plot.data2 <- plyr::join(matm, ave, type = "full", by = "Celltypes")
  scorename <- genesetname
  
  
  if (test.use == "t") {
    stat.test <- matm %>% group_by(features) %>% 
      rstatix::t_test(expression ~celltypes,
                      p.adjust.method = adjust.method, 
                      ref.group = ref.group) %>%
      add_significance( "p.adj")
  }
  if (test.use == "wilcox") {
    stat.test <- matm %>% group_by(features) %>% 
      rstatix::wilcox_test(expression ~celltypes,
                           p.adjust.method = adjust.method, 
                           ref.group = ref.group) %>%
      add_significance( "p.adj")
  }
  
  if (test.use != "none") {
    group_num <- stat.test %>% group_by(features) %>% summarize(count = n())
    y.position.adjust <- unlist(as.data.frame(sapply(group_num$count, 
                                                     function(x) 1:x)))
    y.position <- y.position.adjust * max(matm$expression)/10 * 
      adjust_p_text_gap + max(matm$expression)
    y.position <- y.position * adjust_p_text_height
    
    if (is.null(ylim)) {
      ylim = c(0, max(y.position) *ylim.adjust)
    }
    p <- ggpubr::ggboxplot(matm, x = "celltypes", y = "expression",
                           color = "celltypes", 
                           facet.by = "features", ylim = ylim) + xlab("") 
    
    p <- p + ggpubr::stat_pvalue_manual(stat.test, label = plabel, 
                                        y.position = y.position,
                                        remove.bracket = remove.bracket,...)
    
  }else{
    p <- ggpubr::ggboxplot(matm, x = "celltypes", y = "expression",
                           color = "celltypes", 
                           facet.by = "features", ylim = ylim) + xlab("") 
  }
  if (is.null(color.use)) {
    return(p + ggsci::scale_color_jco())
  }
  p <- p + scale_color_manual(values = color.use)
  return(p)
}



#' Calculate gene score
#' 
#' @return Data frame
#' 
#' @export
data_score <- function (object, genesetlist, cellpos) 
{
  datalist <- list()
  for (st in names(genesetlist)) {
    geneset <- genesetlist[[st]]
    Mgene_in_data <- geneset %in% rownames(object)
    missgene <- geneset[!Mgene_in_data]
    if (any(!Mgene_in_data)) {
      message(paste0(st, " signatures: ", paste(missgene, 
                                                collapse = ","), " are not in data."))
    }
    mgene <- geneset[Mgene_in_data]
    plot.data <- as.data.frame(t(as.matrix(object@assays$RNA@data[mgene, 
                                                                  cellpos])))
    plot.data <- as.data.frame(apply(plot.data, 1, mean))
    plot.data$Celltypes <- as.vector(Idents(object))[cellpos]
    plot.data$Barcodes <- names(Idents(object))[cellpos]
    colnames(plot.data) <- c("Expression", "Celltypes", "Barcodes")
    ave <- plot.data %>% group_by(Celltypes) %>% summarize(mean_exp = mean(Expression), 
                                                           .groups = "drop")
    plot.data2 <- plyr::join(plot.data, ave, type = "full", 
                             by = "Celltypes")
    colnames(plot.data2)[c(1, 4)] <- paste0(st, "_", colnames(plot.data2)[c(1, 
                                                                            4)])
    datalist[[st]] <- plot.data2
  }
  gname <- names(genesetlist)
  plotdata <- inner_join(datalist[[gname[1]]], datalist[[gname[2]]], 
                         by = c("Celltypes", "Barcodes"))
  return(plotdata)
}
