#' Volcano plot
#' 
#' @export
volcanoplot <- function (DEmat, avg_logFC = "avg_logFC", p_val_adj = "p_val_adj", 
    feature = "gene", color.use = c("red", "grey"), p_threshold = 0.05, 
    FC_threshold = 0.25, topn = 5, show_gene = NULL, only.show.sign.text = T,
    text.size = 5, pt.size = 1, title = "") 
{
    options(ggrepel.max.overlaps = 10)
    DEmat <- data.frame(feature = DEmat[, feature], avg_logFC = DEmat[, 
        avg_logFC], p_val_adj = DEmat[, p_val_adj], stringsAsFactors = F)
    DEmat$cluster <- "Pos"
    DEmat$cluster[DEmat$avg_logFC < 0] <- "Neg"
    DEmat$Significant <- ifelse(DEmat$p_val_adj < p_threshold & 
        abs(DEmat$avg_logFC) > FC_threshold, paste0("FDR < ", 
        p_threshold, " and log2FC > ", FC_threshold), "Not Significant")
    
    if (topn > 0) {
        top_p_gene <- DEmat %>% group_by(cluster) %>% top_n(topn, 
            -abs(p_val_adj))
        DEmat$show_gene_pos <- DEmat$feature %in% top_p_gene$feature
    }else{
      DEmat$show_gene_pos <- FALSE
    }
    if (!is.null(show_gene)) {
        DEmat$show_gene_pos <- DEmat$show_gene_pos | (DEmat$feature %in% 
            show_gene)
    }
    if (only.show.sign.text) {
      DEmat$show_gene_pos <- DEmat$Significant != "Not Significant" & 
            DEmat$show_gene_pos
    }
   # 
    DEmat$feature[!DEmat$show_gene_pos] <- ""
    if (sum(DEmat$Significant == "Not Significant") == nrow(DEmat)) {
        color.use <- "grey"
    }
    p <- ggplot(DEmat, aes(x = avg_logFC, y = -log10(p_val_adj))) + 
        geom_point(aes(color = Significant), size = pt.size) + 
        scale_color_manual(values = color.use) + geom_vline(aes(xintercept = FC_threshold), 
        colour = "black", linetype = "dashed") + geom_vline(aes(xintercept = -FC_threshold), 
        colour = "black", linetype = "dashed") + theme_bw() + 
        xlab("log2(Fold Change)") + ylab("-log10(FDR)") + ggtitle(title) + 
        theme(legend.position = "bottom", panel.grid.major = element_line(colour = NA), 
            panel.grid.minor = element_line(colour = NA), plot.title = element_text(size = rel(1.5), 
                hjust = 0.5)) + ggrepel::geom_text_repel(aes(x = avg_logFC, 
        y = -log10(p_val_adj), label = feature), size = text.size, 
        max.overlaps = Inf)
    p
    return(p)
}


#' Diagonal volcano plot
#' 
#' @export
volcanoplot2 <- function (object, DEmat, annotation, ident.1, ident.2, assay = "RNA", 
    avg_logFC_name = "avg_logFC", p_val_adj_name = "p_val_adj",
    legend_logFC_name = 'logFC',
    legend_p_value_name = 'FDR',
    feature_name = "gene", do.expm1 = F, 
    color.use = c("red", "grey"),
    p_threshold = 0.05, FC_threshold = 0.25, 
    show.gene = F, topn_p = NULL, topn_FC = NULL, add.gene.list = NULL, 
    AI.friendly = F,
    text.size = 3, pt.size = 1, legend.point.size = 3, 
    repel.max.overlaps = Inf, min.segment.length = 0.5, max.iter = 10000,
    title = "",...) {
  
    options(ggrepel.max.overlaps = 10)
    DEmat <- data.frame(gene = as.vector(DEmat[, feature_name]), 
                        avg_logFC = DEmat[, avg_logFC_name],
                        p_val_adj = DEmat[, p_val_adj_name], 
                        stringsAsFactors = F)
    DEmat$Significant <- ifelse(DEmat$p_val_adj < p_threshold & 
        (abs(DEmat$avg_logFC) > FC_threshold), paste0(legend_p_value_name," < ", 
        p_threshold, " and ",legend_logFC_name," > ", FC_threshold),
        "Not Significant")
    DEmat$cluster <- paste(ident.1, collapse = "-")
    DEmat$cluster[DEmat$avg_logFC < 0] <- paste(ident.2, collapse = "-")
    if (show.gene) {
      DEmat$show_gene_pos <- F
      if (!is.null(topn_p)) {
        topn_gene <- DEmat %>% top_n(topn_p, -abs(p_val_adj))
        DEmat$show_gene_pos <- DEmat$gene %in% topn_gene$gene
      }
      if (!is.null(topn_FC)) {
        topn_gene <- DEmat %>% group_by(cluster) %>% top_n(topn_FC, 
                                                           abs(avg_logFC))
        DEmat$show_gene_pos <- DEmat$gene %in% topn_gene$gene
      }
      
      if (!is.null(add.gene.list)) {
        DEmat$show_gene_pos <- DEmat$show_gene_pos | (DEmat$gene %in% 
                                                        add.gene.list)
      }
      DEmat$label <- DEmat$gene
      DEmat$label[!DEmat$show_gene_pos] <- ""
    }else{
      DEmat$label <- ''
    }
    df1 <- object@assays[[assay]]@data[DEmat$gene, object@meta.data[, 
        annotation] %in% ident.1]
    df2 <- object@assays[[assay]]@data[DEmat$gene, object@meta.data[, 
        annotation] %in% ident.2]
    if (do.expm1) {
        df1 <- expm1(df1)
        df2 <- expm1(df2)
    }
    plotdata <- data.frame(x = Matrix::rowMeans(df1), y = Matrix::rowMeans(df2), 
        Significant = DEmat$Significant, label = DEmat$label, 
        stringsAsFactors = F)
    p <- ggplot(plotdata, aes(x, y, color = Significant)) + 
      geom_point(size = pt.size) + 
      scale_color_manual(values = color.use) + 
      theme_bw() + 
      xlab(paste0("Mean expression of ", paste(ident.1, collapse = "-"))) + 
      ylab(paste0("Mean expression of ", paste(ident.2, collapse = "-"))) + 
      ggtitle(title) + theme(legend.position = "bottom", 
                             panel.grid.major = element_line(color = NA),
                             panel.grid.minor = element_line(color = NA), 
                             plot.title = element_text(size = rel(1.5), 
                                                       hjust = 0.5))
    if (AI.friendly) {
      p <- ggAIplot(p)
    }
    p <- p + ggrepel::geom_text_repel(aes(x, y, label = label), 
                                      size = text.size,
                                      color = 'black',
                                      max.iter = max.iter,
                                      min.segment.length = min.segment.length,
                                      max.overlaps = repel.max.overlaps,
                                      ...)
    p <- p + guides(color = guide_legend(override.aes = list(size=legend.point.size)))
    return(p)
}
