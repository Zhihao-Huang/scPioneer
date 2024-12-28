#' Dotplot of GO enrichment.
#'
#' @param GOmat Result of GO enrichment from clusterGO.
#' @param show.top.term Number of top GO terms of each cluster showed in dotplot.
#' @param log.adj.P Show -log(adjust p value) instead of adjust p value.
#' @param cell.order Cell order showed in dotplot.
#' @param term.order GO terms order showed in dotplot.
#' @param coord_flip Parameter in ggplot.
#' @param text.y.size Text size of GO terms.
#' @param text.x.size Text size of Celltypes.
#' @param color.use Color patterns used for showing p value.
#' @param legend.title.size Text size of legend title.
#' @param legend.text.size Text size of legend text.
#' @examples
#' markerlist <- FindAllMarkers(pbmc,only.pos = TRUE,logfc.threshold = 0.5)
#' GOmat <- clusterGO(markerlist)
#' GOplot(GOmat)
#' @export
GOplot <- function(
  GOmat,
  show.top.term = 5,
  log.adj.P = F,
  cell.order = NULL,
  term.order = NULL,
  coord_flip = F,
  color.use = NULL,
  text.y.size = 8,
  text.x.size = 10,
  legend.title.size = 10,
  legend.text.size = 10){
  topnGO <- GOmat %>% group_by(celltype) %>% top_n(show.top.term,-p.adjust)
  plot.go.id<- unique(topnGO$ID)
  plotdata <- GOmat[GOmat$ID %in% plot.go.id,]
  plotdata <- plotdata[,c('Description','p.adjust','celltype',"GeneRatio")]
  plotdata$log.p.adjust <- round(-log(plotdata$p.adjust),2)
  plotdata$log.p.adjust[plotdata$p.adjust > 0.05] <- 0
  plotdata$GeneRatio <- round(as.numeric(gsub('\\/.*$','',
                                              plotdata$GeneRatio)) / as.numeric(gsub('^.*\\/','',plotdata$GeneRatio)),2)
  colnames(plotdata) <- c('GO_term','adjust.P','Celltype','GeneRatio','log.adj.P')
  if(!is.null(cell.order)){
    plotdata$Celltype <- factor(plotdata$Celltype,levels = cell.order)
  }else{
    plotdata$Celltype <- factor(plotdata$Celltype,levels = unique(topnGO$celltype))
  }
  if(!is.null(term.order)){
    plotdata$GO_term <- factor(plotdata$GO_term,levels = term.order)
  }else{
    plotdata$GO_term <- factor(plotdata$GO_term,levels = rev(unique(topnGO$Description)))
  }
  if (!is.null(color.use)){
    color.use <- color.use
  }else{
    color.use <- RColorBrewer::brewer.pal(11,"RdYlBu")
  }
  if (log.adj.P) {
    p <- ggplot(plotdata, aes(x=Celltype, y=GO_term, colour = log.adj.P, size=GeneRatio))
    color.use <- rev(RColorBrewer::brewer.pal(11,"RdYlBu"))
  }else{
    p <- ggplot(plotdata, aes(x=Celltype, y=GO_term, colour = adjust.P, size=GeneRatio))
  }
  mytheme <- theme(#plot.title = element_blank(),
                   #element_text(family="ArialMT", size=25, hjust = 0.5),
                   #panel.grid= element_blank(),
                   axis.title = element_blank(),
                   axis.ticks.length = unit(0.3, "cm"),
                   axis.text=element_text(family="ArialMT", size = text.y.size),
                   axis.text.x=element_text(angle=90, vjust=0.5, hjust=1,size= text.x.size),
                   legend.title=element_text(family="ArialMT", size=legend.title.size),
                   legend.text=element_text(family="ArialMT", size=legend.text.size),
                   legend.background=element_rect(fill="transparent"),
                   plot.margin = unit(c(1,1,1,5), "cm"))
  p <- p + geom_point() +
    scale_color_gradientn(name='adjust P value', colors=color.use) +
    guides(size=guide_legend(title="GeneRatio")) +theme_bw() +
    mytheme
  if (coord_flip) {
    p <- p + coord_flip()
  }
  return(p)
}

#' Heatmap of GO enrichment.
#'
#' @param GOmat Result of GO enrichment from clusterGO.
#' @param show.top.term Number of top GO terms of each cluster showed in dotplot.
#' @param log.adj.P Show -log(adjust p value) instead of adjust p value.
#' @param cell.order Cell order showed in dotplot.
#' @param term.order GO terms order showed in dotplot.
#' @param coord_flip Parameter in ggplot.
#' @param text.y.size Text size of GO terms.
#' @param text.x.size Text size of Celltypes.
#' @param color.use Color patterns used for showing p value.
#' @param legend.title.size Text size of legend title.
#' @param legend.text.size Text size of legend text.
#' @examples
#' markerlist <- FindAllMarkers(pbmc,only.pos = TRUE,logfc.threshold = 0.5)
#' GOmat <- clusterGO(markerlist)
#' GOheatmap(GOmat)
#' @export
GOheatmap <- function(
  GOmat,
  show.top.term = 5,
  log.adj.P = F,
  cell.order = NULL,
  term.order = NULL,
  coord_flip = F,
  color.use = NULL,
  text.y.size = 8,
  text.x.size = 10,
  legend.title.size = 10,
  legend.text.size = 10,
  method = c('pheatmap','ggplot2')[2],
  cluster_rows=FALSE,
  cluster_cols=FALSE){
  if (is.null(GOmat$ID)) {
    GOmat$ID <- GOmat$Description
  }
  if (is.null(GOmat$celltype) & 'cluster' %in% colnames(GOmat) ) {
    GOmat$celltype <- GOmat$cluster
  }
  if (is.null(GOmat$celltype) & 'Celltype' %in% colnames(GOmat) ) {
    GOmat$celltype <- GOmat$Celltype
  }
  topnGO <- GOmat %>% group_by(celltype) %>% top_n(show.top.term,-p.adjust)
  plot.go.id<- unique(topnGO$ID)
  plotdata <- topnGO[topnGO$ID %in% plot.go.id,]
  plotdata <- plotdata[,c('Description','p.adjust','celltype',"GeneRatio")]
  plotdata$log.p.adjust <- round(-log(plotdata$p.adjust),2)
  plotdata$log.p.adjust[plotdata$p.adjust > 0.05] <- 0
  plotdata$GeneRatio <- round(as.numeric(gsub('\\/.*$','',
                                              plotdata$GeneRatio)) / as.numeric(gsub('^.*\\/','',plotdata$GeneRatio)),2)
  colnames(plotdata) <- c('GO_term','adjust.P','Celltype','GeneRatio','log.adj.P')
  plotdata1 <- plotdata[,c('GO_term','Celltype','adjust.P')]
  pos <- duplicated(plotdata1[,c('GO_term','Celltype')])
  if(sum(pos) > 0) {
    message('Duplicated pathways in same celltype were found.')
  }
  plotdata1 <- plotdata1[!pos,]
  ddata <- dcast(plotdata1, GO_term~Celltype,value.var = 'adjust.P')
  ddata[is.na(ddata)] <- 0.05
  matm <- melt(ddata, id.vars = 'GO_term', measure.vars = colnames(ddata)[-1], variable.name = 'Celltype',value.name = 'adjust.P')
  rownames(ddata) <- ddata$GO_term
  ddata <- ddata[-1]
  if(!is.null(cell.order)){
    pos <- cell.order %in% colnames(ddata)
    if (sum(!pos) > 0) {
      warning('Cell types from cell_order are not in enrichment data: \n', paste0(cell.order[!pos],collapse = ', '))
    }
    cell.order <-  cell.order[pos]
    plotdata$Celltype <- factor(plotdata$Celltype,levels = cell.order)
    matm$Celltype <- factor(matm$Celltype,levels = cell.order)
    matm <- matm %>% arrange(Celltype, adjust.P)
    topnGO$celltype <- factor(topnGO$celltype, levels = cell.order)
    ddata <- ddata[,cell.order]
  }else{
    plotdata$Celltype <- factor(plotdata$Celltype,levels = unique(topnGO$celltype))
    matm$Celltype <- factor(matm$Celltype,levels = unique(topnGO$celltype))
    ddata <- ddata[,unique(as.vector(topnGO$celltype))]
  }
  if(!is.null(term.order)){
    pos <- term.order %in% rownames(ddata)
    if (sum(!pos) > 0) {
      warning('Terms from term.order are not in enrichment data: \n', paste0(term.order[!pos],collapse = ', '))
    }
    term.order <- term.order[pos]
    plotdata$GO_term <- factor(plotdata$GO_term,levels = term.order)
    matm$GO_term <- factor(matm$GO_term,levels = term.order)
    ddata <- ddata[term.order,]
  }else{
    plotdata$GO_term <- factor(plotdata$GO_term,levels = rev(unique(topnGO$Description)))
    topnGO <- topnGO %>% arrange(celltype, p.adjust)
    matm$GO_term <- factor(matm$GO_term,levels = rev(unique(topnGO$Description)))
    ddata <- ddata[rev(unique(as.vector(topnGO$Description))),]
  }
  if (!is.null(color.use)){
    color.use1 <- color.use
    color.use2 <- color.use
  }else{
    color.use2 <- c(RColorBrewer::brewer.pal(11,"RdYlBu")[1], RColorBrewer::brewer.pal(11,"RdYlBu")[3],RColorBrewer::brewer.pal(11,"RdYlBu")[6:11])
    color.use1 <- c(RColorBrewer::brewer.pal(11,"RdYlBu")[2],"white","steelblue")
    
  }
  if (method == 'ggplot2') {
    mytheme <- theme(plot.title = element_blank(),
                     #element_text(family="ArialMT", size=25, hjust = 0.5),
                     #panel.grid= element_blank(),
                     axis.title = element_blank(),
                     axis.ticks.length = unit(0.3, "cm"),
                     axis.text=element_text(family="ArialMT", size = text.y.size),
                     axis.text.x=element_text(angle=90, vjust=0.5, hjust=1,size= text.x.size),
                     legend.title=element_text(family="ArialMT", size=legend.title.size),
                     legend.text=element_text(family="ArialMT", size=legend.text.size),
                     legend.background=element_rect(fill="transparent"),
                     plot.margin = unit(c(1,1,1,5), "cm"))
    if (log.adj.P) {
      p <- ggplot(matm, aes(Celltype, GO_term)) +
        geom_tile(aes(fill = log.adj.P),colour = "black")+ scale_fill_gradientn(colors = color.use2)+
        mytheme
    }else{
      p <- ggplot(matm, aes(Celltype, GO_term)) +
        geom_tile(aes(fill = adjust.P),colour = "black")+ scale_fill_gradientn(colors = color.use2)+
        mytheme
    }
    return(p)
  }else if (method == 'pheatmap') {
    return(pheatmap::pheatmap(ddata, color = color.use2,
                              cluster_rows=cluster_rows, cluster_cols=cluster_cols))
  }
}

#' bar plot for GO enrichemnt of single group.
#' 
#' @export
GObarplot <- function (data) {
  data$terms <- paste0(data$ID, ": ", data$Description)
  data$logadjp <- -log10(data$p.adjust)
  data$logp <- -log10(data$pvalue)
  data <- data[order(data$logp), ]
  data$terms <- factor(data$terms, levels = unique(data$terms))
  p <- ggplot(data, aes(terms, logp, fill = logp)) + geom_hline(yintercept = c(2, 
                                                                               4), color = "grey") + geom_bar(stat = "identity", width = 0.7) + 
    coord_flip() + theme_bw() + theme(panel.grid.major = element_line(colour = NA), 
                                      panel.grid.minor = element_line(colour = NA), axis.ticks.y = element_blank(), 
                                      legend.position = "none") + scale_x_discrete(position = "top") + 
    scale_fill_gradientn(colours = c("#F4A460", "#D2691E", 
                                     "#8B4513")) + ylab("-log10(P)") + xlab("")
  return(p)
}

#' Heatmap of GO enrichment.
#'
#' @param GOmat Result of GO enrichment from clusterGO.
#' @param show.top.term Number of top GO terms of each cluster showed in dotplot.
#' @param log.adj.P Show -log(adjust p value) instead of adjust p value.
#' @param cell.order Cell order showed in dotplot.
#' @param term.order GO terms order showed in dotplot.
#' @param coord_flip Parameter in ggplot.
#' @param text.y.size Text size of GO terms.
#' @param text.x.size Text size of Celltypes.
#' @param color.use Color patterns used for showing p value.
#' @param legend.title.size Text size of legend title.
#' @param legend.text.size Text size of legend text.
#' @examples
#' markerlist <- FindAllMarkers(pbmc,only.pos = TRUE,logfc.threshold = 0.5)
#' GOmat <- clusterGO(markerlist)
#' GOheatmap(GOmat)
#' @export
GOheatmap <- function (GOmat, show.top.term = 5, log10.adj.P = F, na.as = 0.05, 
    cell.order = NULL, term.order = NULL, coord_flip = F, color.use = NULL, 
    grid.color = "black", text.y.size = 8, text.x.size = 10, 
    legend.title.size = 10, legend.text.size = 10, max.length = 35, 
    method = c("pheatmap", "ggplot2")[2], cluster_rows = FALSE, 
    cluster_cols = FALSE) 
{
    if (is.null(GOmat$ID)) {
        GOmat$ID <- GOmat$Description
    }
    if (is.null(GOmat$celltype) & "cluster" %in% colnames(GOmat)) {
        GOmat$celltype <- GOmat$cluster
    }
    if (is.null(GOmat$celltype) & "Celltype" %in% colnames(GOmat)) {
        GOmat$celltype <- GOmat$Celltype
    }
    topnGO <- GOmat %>% group_by(celltype) %>% top_n(show.top.term, 
        -p.adjust)
    plot.go.id <- unique(topnGO$ID)
    plotdata <- topnGO[topnGO$ID %in% plot.go.id, ]
    plotdata <- plotdata[, c("Description", "p.adjust", "celltype", 
        "GeneRatio")]
    plotdata$log.p.adjust <- round(-log10(plotdata$p.adjust), 2)
    plotdata$log.p.adjust[plotdata$p.adjust > na.as] <- 0
    plotdata$GeneRatio <- round(as.numeric(gsub("\\/.*$", "", 
        plotdata$GeneRatio))/as.numeric(gsub("^.*\\/", "", plotdata$GeneRatio)), 
        2)
    colnames(plotdata) <- c("GO_term", "adjust.P", "Celltype", 
        "GeneRatio", "log.adj.P")
    plotdata1 <- plotdata[, c("GO_term", "Celltype", "adjust.P")]
    pos <- duplicated(plotdata1[, c("GO_term", "Celltype")])
    if (sum(pos) > 0) {
        message("Duplicated pathways in same celltype were found.")
    }
    plotdata1 <- plotdata1[!pos, ]
    ddata <- dcast(plotdata1, GO_term ~ Celltype, value.var = "adjust.P")
    ddata[is.na(ddata)] <- na.as
    matm <- melt(ddata, id.vars = "GO_term", measure.vars = colnames(ddata)[-1], 
        variable.name = "Celltype", value.name = "adjust.P")
    rownames(ddata) <- ddata$GO_term
    ddata <- ddata[-1]
    if (!is.null(cell.order)) {
        pos <- cell.order %in% colnames(ddata)
        if (sum(!pos) > 0) {
            warning("Cell types from cell_order are not in enrichment data: \n", 
                paste0(cell.order[!pos], collapse = ", "))
        }
        cell.order <- cell.order[pos]
        plotdata$Celltype <- factor(plotdata$Celltype, levels = cell.order)
        matm$Celltype <- factor(matm$Celltype, levels = cell.order)
        matm <- matm %>% arrange(Celltype, adjust.P)
        topnGO$celltype <- factor(topnGO$celltype, levels = cell.order)
        ddata <- ddata[, cell.order]
    }
    else {
        plotdata$Celltype <- factor(plotdata$Celltype, levels = unique(topnGO$celltype))
        matm$Celltype <- factor(matm$Celltype, levels = unique(topnGO$celltype))
        ddata <- ddata[, unique(as.vector(topnGO$celltype))]
    }
    if (!is.null(term.order)) {
        pos <- term.order %in% rownames(ddata)
        if (sum(!pos) > 0) {
            warning("Terms from term.order are not in enrichment data: \n", 
                paste0(term.order[!pos], collapse = ", "))
        }
        term.order <- term.order[pos]
        plotdata$GO_term <- factor(plotdata$GO_term, levels = term.order)
        matm$GO_term <- factor(matm$GO_term, levels = term.order)
        ddata <- ddata[term.order, ]
    }
    else {
        plotdata$GO_term <- factor(plotdata$GO_term, levels = rev(unique(topnGO$Description)))
        topnGO <- topnGO %>% arrange(celltype, p.adjust)
        matm$GO_term <- factor(matm$GO_term, levels = rev(unique(topnGO$Description)))
        ddata <- ddata[rev(unique(as.vector(topnGO$Description))), 
            ]
    }
    if (!is.null(color.use)) {
        color.use1 <- color.use
        color.use2 <- color.use
    }
    else {
        color.use2 <- c(RColorBrewer::brewer.pal(11, "RdYlBu")[1], 
            RColorBrewer::brewer.pal(11, "RdYlBu")[3], RColorBrewer::brewer.pal(11, 
                "RdYlBu")[6:11])
        color.use1 <- rev(scPalette_heatmap_bar())
    }
    if (method == "ggplot2") {
        mytheme <- theme(plot.title = element_blank(), axis.title = element_blank(), 
            axis.ticks.length = unit(0.3, "cm"), axis.text = element_text(family = "ArialMT", 
                size = text.y.size), axis.text.x = element_text(angle = 90, 
                vjust = 0.5, hjust = 1, size = text.x.size), 
            legend.title = element_text(family = "ArialMT", size = legend.title.size), 
            legend.text = element_text(family = "ArialMT", size = legend.text.size), 
            legend.background = element_rect(fill = "transparent"), 
            plot.margin = unit(c(1, 1, 1, 5), "cm"))
        if (log10.adj.P) {
            matm$log.adj.P <- round(-log10(matm$adjust.P), 2)
            matm$log.adj.P[matm$adjust.P > na.as] <- 0
            p <- ggplot(matm, aes(Celltype, GO_term)) + geom_tile(aes(fill = log.adj.P), 
                colour = grid.color) + scale_fill_gradientn(name = "-log10(adjust P)", 
                colors = color.use1) + mytheme
        }
        else {
            p <- ggplot(matm, aes(Celltype, GO_term)) + geom_tile(aes(fill = adjust.P), 
                colour = grid.color) + scale_fill_gradientn(colors = color.use2) + 
                mytheme
        }
        p <- p + scale_y_discrete(labels = function(x) stringr::str_wrap(x, 
            width = max.length))
        return(p)
    }
    else if (method == "pheatmap") {
        return(pheatmap::pheatmap(ddata, color = color.use2, 
            cluster_rows = cluster_rows, cluster_cols = cluster_cols))
    }
}

#' Dotplot of GO enrichment.
#'
#' @param GOmat Result of GO enrichment from clusterGO.
#' @param show.top.term Number of top GO terms of each cluster showed in dotplot.
#' @param log.adj.P Show -log(adjust p value) instead of adjust p value.
#' @param cell.order Cell order showed in dotplot.
#' @param term.order GO terms order showed in dotplot.
#' @param coord_flip Parameter in ggplot.
#' @param text.y.size Text size of GO terms.
#' @param text.x.size Text size of Celltypes.
#' @param color.use Color patterns used for showing p value.
#' @param legend.title.size Text size of legend title.
#' @param legend.text.size Text size of legend text.
#' @examples
#' markerlist <- FindAllMarkers(pbmc,only.pos = TRUE,logfc.threshold = 0.5)
#' GOmat <- clusterGO(markerlist)
#' GOplot(GOmat)
#' @export
GOplot <- function(
  GOmat,
  show.top.term = 5,
  log.adj.P = F,
  cell.order = NULL,
  term.order = NULL,
  coord_flip = F,
  color.use = NULL,
  text.y.size = 8,
  text.x.size = 10,
  legend.title.size = 10,
  legend.text.size = 10){
  topnGO <- GOmat %>% group_by(celltype) %>% top_n(show.top.term,-p.adjust)
  plot.go.id<- unique(topnGO$ID)
  plotdata <- GOmat[GOmat$ID %in% plot.go.id,]
  plotdata <- plotdata[,c('Description','p.adjust','celltype',"GeneRatio")]
  plotdata$log.p.adjust <- round(-log(plotdata$p.adjust),2)
  plotdata$log.p.adjust[plotdata$p.adjust > 0.05] <- 0
  plotdata$GeneRatio <- round(as.numeric(gsub('\\/.*$','',
                                              plotdata$GeneRatio)) / as.numeric(gsub('^.*\\/','',plotdata$GeneRatio)),2)
  colnames(plotdata) <- c('GO_term','adjust.P','Celltype','GeneRatio','log.adj.P')
  if(!is.null(cell.order)){
    plotdata$Celltype <- factor(plotdata$Celltype,levels = cell.order)
  }else{
    plotdata$Celltype <- factor(plotdata$Celltype,levels = unique(topnGO$celltype))
  }
  if(!is.null(term.order)){
    plotdata$GO_term <- factor(plotdata$GO_term,levels = term.order)
  }else{
    plotdata$GO_term <- factor(plotdata$GO_term,levels = rev(unique(topnGO$Description)))
  }
  if (!is.null(color.use)){
    color.use <- color.use
  }else{
    color.use <- RColorBrewer::brewer.pal(11,"RdYlBu")
  }
  if (log.adj.P) {
    p <- ggplot(plotdata, aes(x=Celltype, y=GO_term, colour = log.adj.P, size=GeneRatio))
    color.use <- rev(RColorBrewer::brewer.pal(11,"RdYlBu"))
  }else{
    p <- ggplot(plotdata, aes(x=Celltype, y=GO_term, colour = adjust.P, size=GeneRatio))
  }
  mytheme <- theme(plot.title = element_blank(),
                   #element_text(family="ArialMT", size=25, hjust = 0.5),
                   #panel.grid= element_blank(),
                   axis.title = element_blank(),
                   axis.ticks.length = unit(0.3, "cm"),
                   axis.text=element_text(family="ArialMT", size = text.y.size),
                   axis.text.x=element_text(angle=90, vjust=0.5, hjust=1,size= text.x.size),
                   legend.title=element_text(family="ArialMT", size=legend.title.size),
                   legend.text=element_text(family="ArialMT", size=legend.text.size),
                   legend.background=element_rect(fill="transparent"),
                   plot.margin = unit(c(1,1,1,5), "cm"))
  p <- p + geom_point() +
    scale_color_gradientn(name='adjust P value', colors=color.use) +
    guides(size=guide_legend(title="GeneRatio")) +theme_bw() +
    mytheme
  if (coord_flip) {
    p <- p + coord_flip()
  }
  return(p)
}

#' Bar plot to show the result of all clusters' GESA.
#' 
#' @examples
#' DEG_GSEA <- parallel::mclapply(levels(Idents(pbmc)), function(x) {
#' mmk <- FindMarkers(pbmc, ident.1 = x, test.use = 'wilcox',
#'                    min.pct = 0.25,logfc.threshold = 0)
#' mmk$gene <- rownames(mmk)
#' mmk$cluster <- x
#' mmk$cluster[mmk$avg_logFC < 0] <- paste0(x,'_Others')
#' mmk
#' }, mc.cores = 10)
#' names(DEG_GSEA) <- levels(Idents(pbmc))
#' gsealist <- GSEA_all(DEG_GSEA, category = 'C2',subcategory = NULL, n.cores = 10)
#' fgsea_df <- do.call(rbind, gsealist)
#' GSEA_all_plot(fgsea_df, show.num = 3)
#' 
#' @export
GSEA_all_plot <- function (fgsea_df, only.pos = T, color.use = NULL, show.num = 15, 
    pvalue = 0.05, text.gap = 0.1, title = NULL) 
{
    if (only.pos) {
        fgsea_df <- fgsea_df[fgsea_df$NES > 0, ]
    }
    df <- fgsea_df %>% group_by(celltype) %>% filter(pval < pvalue) %>% 
        arrange(desc(abs(NES))) %>% top_n(n = show.num, wt = NES) %>% 
        mutate(pathway = tidytext::reorder_within(pathway, NES, 
            celltype))
    df$hjust <- ifelse(df$NES > 0, 1, 0)
    df$ypos <- ifelse(df$NES < 0, text.gap, -text.gap)
    if (is.null(color.use)) {
        color.use <- scPalette2(length(unique(df$celltype)))
    }
    p <- ggplot(df, aes(pathway, NES, label = pathway, hjust = hjust)) + 
        geom_col(aes(fill = celltype)) + coord_flip() + facet_grid(rows = vars(celltype), 
        scales = "free", space = "free") + labs(x = "", y = "Normalized Enrichment Score", 
        title = title) + theme_classic() + scale_fill_manual(values = color.use, 
        name = "Cell types") + theme(legend.position = "top", 
        strip.text = element_blank())
    p <- p + tidytext::scale_x_reordered()
    if (any(df$NES < 0)) {
        p <- p + geom_text(aes(y = ypos), colour = "black", size = 3) + 
            scale_x_discrete(breaks = NULL) + theme(axis.line.y = element_blank())
    }
    return(p)
}

#' pairs barplot showing the result of GESA.
#' 
#' @examples 
#' celltype_1 <- 'CD8 T'
#' celltype_2 <- 'NK'
#' DEG_pair <- FindMarkers(pbmc, ident.1 = celltype_1, ident.2 = celltype_2, min.pct = 0.25, logfc.threshold = 0)
#' DEG_pair$gene <- rownames(DEG_pair)
#' DEG_pair$cluster <- celltype_1
#' DEG_pair$cluster[DEG_pair$avg_logFC < 0] <- celltype_2
#' 
#' GSEA_KEGG <- GSEA_pair(DEG_pair, category = 'C2', subcategory = 'KEGG')
#' GSEA_KEGG <- GSEA_KEGG[GSEA_KEGG$padj < 0.05, ]
#' GSEA_C2 <- GSEA_pair(DEG_pair, category = 'C2', subcategory = NULL)
#' GSEA_C2 <- GSEA_C2[GSEA_C2$padj < 0.05, ]
#' 
#' GSEA_pair_plot(GSEA_C2)
#' 
#' @export
GSEA_pair_plot <- function (fgsea_df, show.num = 15, pvalue.adj = 0.05,
                            text.gap = 0.1, text.max.length = 50, 
                            bar.width = 0.8, text.lineheight = 0.6, text.size = 3,
                            show.gene.on.bar = F, gene.max.length = 10, 
                            gene.size = 3, gene.gap = 0.1,
                            color = NULL, 
                            title = NULL, legend.title = NULL) {
  df <- fgsea_df %>% group_by(celltype) %>% arrange(desc(abs(NES))) %>% 
    filter(padj < pvalue.adj) %>% top_n(n = show.num, wt = NES)
  df$hjust <- ifelse(df$NES > 0, 1, 0)
  df$hjust.gene <- ifelse(df$NES > 0, 0, 1)
  df$ypos <- ifelse(df$NES < 0, text.gap, -text.gap)
  df$ypos.gene <- ifelse(df$NES < 0, -gene.gap, gene.gap)
  if (is.null(color)) {
    color <- scPalette2(2)
  }
  df$pathway <- sapply(df$pathway, function(x) text_multi_lines(x,text.max.length = text.max.length))
  p <- ggplot(df, aes(reorder(pathway, NES), NES, label = pathway, 
                      hjust = hjust)) + 
    geom_col(aes(fill = celltype), width = bar.width) + coord_flip() + 
    labs(x = "", y = "Normalized Enrichment Score", title = title) + 
    theme_classic() + scale_fill_manual(values = color, 
                                        name = legend.title) + 
    theme(legend.position = "top")
  if (any(df$NES < 0)) {
    p <- p + geom_text(aes(y = ypos), 
                       colour = "black", size = text.size, 
                       lineheight = text.lineheight) + 
      scale_x_discrete(breaks = NULL) + theme(axis.line.y = element_blank())
  }
  if (show.gene.on.bar) {
    p$data$leadingEdge <- sapply(p$data$leadingEdge, function(x) {
      text_ellipsis(x,gene.max.length = gene.max.length)})
    p <- p + geom_text(aes(y = ypos.gene, label = leadingEdge, hjust = hjust.gene), 
                       colour = "black", size = gene.size, 
                       lineheight = text.lineheight)
      
  }
  return(p)
}

#' Convert string to multi lines.
#' 
#' @export
text_multi_lines <- function(string, text.max.length = 50) {
  char_num <- nchar(string)
  if (char_num < text.max.length) {
    return(string)
  }
  strlist <- strsplit(string, split = ' ')[[1]]
  len_per_word <- sapply(strlist, function(x) nchar(x))
  new_string <- c()
  cumlen <- 0
  for(i in 1:length(strlist)) {
    cumlen <- cumlen + len_per_word[i]
    if (cumlen < text.max.length) {
      new_string <- c(new_string, strlist[i])
      new_string <- paste(new_string, collapse = ' ')
    }else{
      new_string <- paste0(new_string, '\n', strlist[i])
      cumlen <- len_per_word[i]
    }
  }
  return(new_string)
}
#' Convert string to ellipsis.
#' 
#' @export
text_ellipsis <- function(string, gene.max.length = 5, sep = ',', replaceby = '...') {
  genelist <- strsplit(string, split = sep)[[1]]
  gene_num <- length(genelist)
  if (gene_num > gene.max.length) {
    genelist <- c(genelist[1:gene.max.length], replaceby)
  }
  new_string <- paste(genelist, collapse = sep)
  return(new_string)
}

#' GO bar plot
#' 
#' @export
GO_all_plot <- function(GOmat,color.use = NULL,
                        show.num = 15, p.adjust.cutoff = 0.05, text.gap = 0.1,
                        top.p.adjust = 1e-30,
                        title = NULL
) {
  GOmat$p.adjust[GOmat$pvalue < top.p.adjust] <- top.p.adjust
  GOmat$log10_p_adj <- -log10(GOmat$p.adjust)
  df <- GOmat %>% group_by(celltype) %>% filter(p.adjust < p.adjust.cutoff) %>%
    arrange(desc(log10_p_adj)) %>% top_n(n = show.num, wt = log10_p_adj) %>%
    mutate(Description = tidytext::reorder_within(Description, log10_p_adj, celltype))
  df$hjust <- ifelse(df$log10_p_adj > 0, 1, 0)
  df$ypos <- ifelse(df$log10_p_adj < 0, text.gap, -text.gap)
  if (is.null(color.use)) {
    color.use <- scPalette2(length(unique(df$celltype)))
  }
  p <- ggplot(df, aes(Description,log10_p_adj, label = Description, hjust = hjust)) +
    geom_col(aes(fill = celltype)) +
    coord_flip() +
    facet_grid(rows = vars(celltype), scales = "free",space = "free") +
    labs(x="", y="-log10(adjust P value)",title = title) +
    theme_classic() +
    scale_fill_manual(values = color.use, name = 'Cell types') +
    theme(legend.position = "top",
          strip.text = element_blank())
  p <- p + tidytext::scale_x_reordered()
  if (any(df$log10_p_adj < 0)) {
    p <- p + geom_text(aes(y = ypos), colour = 'black', size = 3) +
      scale_x_discrete(breaks = NULL) +
      theme(axis.line.y = element_blank())
  }
  return(p)
}
#' pairs barplot showing the result of GO.
#' 
#' 
#' @export
GO_pair_plot <- function (GOmat,
                          show.num = 15, pvalue.adj = 0.05, text.gap = 0.1,
                          top.p.adjust = 1e-30, text.max.length = 50, 
                          bar.width = 0.8, text.lineheight = 0.6, text.size = 3,
                          show.gene.on.bar = F, gene.max.length = 10, 
                          gene.size = 3, gene.gap = 0.1,
                          color = NULL, 
                          plot.title = NULL, legend.title = NULL) {
  GOmat$p.adjust[GOmat$pvalue < top.p.adjust] <- top.p.adjust
  GOmat$log10_p_adj <- -log10(GOmat$p.adjust)
  pos <- GOmat$celltype == unique(GOmat$celltype)[1]
  GOmat$log10_p_adj[pos] <- -GOmat$log10_p_adj[pos]
  df <- GOmat %>% group_by(celltype) %>% 
    filter(p.adjust < pvalue.adj) %>%
    top_n(n = show.num, wt = abs(log10_p_adj)) %>% 
    arrange(desc(abs(log10_p_adj)))
  df <- df[order(df$celltype, df$p.adjust),]
  df$Description <- as.vector(df$Description)
  df$hjust <- ifelse(df$log10_p_adj > 0, 1, 0)
  df$ypos <- ifelse(df$log10_p_adj < 0, text.gap, -text.gap)
  df$hjust.gene <- ifelse(df$log10_p_adj > 0, 0, 1)
  df$ypos.gene <- ifelse(df$log10_p_adj < 0, -gene.gap, gene.gap)
  if (is.null(color)) {
    color <- scPalette2(2)
  }
  df$Description <- sapply(df$Description, function(x) {
    text_multi_lines(as.character(x), text.max.length = text.max.length)
    })
  df$Description <- factor(df$Description, levels = rev(unique(df$Description)))
  p <- ggplot(df, aes(Description, log10_p_adj, label = Description, 
                      hjust = hjust)) + 
    geom_col(aes(fill = celltype), width = bar.width) + coord_flip() + 
    labs(x = "", y = "-log10(adjust P value)", title = plot.title) + 
    theme_classic() + scale_fill_manual(values = color, 
                                        name = legend.title) + 
    theme(legend.position = "top")
  if (any(df$log10_p_adj < 0)) {
    p <- p + geom_text(aes(y = ypos), 
                       colour = "black", size = text.size, 
                       lineheight = text.lineheight) + 
      scale_x_discrete(breaks = NULL) + theme(axis.line.y = element_blank())
  }
  if (show.gene.on.bar) {
    p$data$geneID <- sapply(p$data$geneID, function(x) {
      text_ellipsis(x,gene.max.length = gene.max.length, sep = '/')})
    p <- p + geom_text(aes(y = ypos.gene, label = geneID, hjust = hjust.gene), 
                       colour = "black", size = gene.size, 
                       lineheight = text.lineheight)
    
  }
  return(p)
}

