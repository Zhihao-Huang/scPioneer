#' Detect similarity of different cell types using variable genes from integration data.
#' 
#' @param var.genes.method Methods for high variable genes selection. var: var.genes from Seurat::FindVariableFeatures; pca: genes with higher SD in PCA dimensions.
#' @param clustering_distance Methods for measuring distance in matrix. cormat: construct cell x cell correlation matrix for pheatmap; pearson: distance defined as (1-Pearson correlation coefficient)/2.
#' @param clustering_distance_cols Methods for measuring distance of cells. Pass parameter to dist() function. binary is also called Jaccard.
#' @param clustering_distance_rows similar to clustering_distance_cols.
#' @param reverse.cor Calculate distance for hierarchical clustering. Refer to https://doi.org/10.1016/j.cell.2021.01.010.
#' @examples 
#' myeloid.integrated <- readRDS('path_to_your_seurat_integrated.rds')
#' datalist <- Find_similar(integrated.object = myeloid.integrated, assay = 'integrated',
#'                         var.genes.method = 'var', Annotation = 'Annotation', plot = F, return.data = T)
#' #correlation heatmap
#' Find_similar(integrated.object = myeloid.integrated, assay = 'integrated',layout = 'heatmap',
#' var.genes.method = 'var', Annotation = 'Annotation', plot = T, return.data = F)
#' #dendrogram
#' Find_similar(integrated.object = myeloid.integrated, assay = 'integrated',layout = 'horizontal.dendrogram',
#' var.genes.method = 'var', Annotation = 'Annotation', k_clusters = 4, plot = T, return.data = F)
#' # cirlize denrogram
#' Find_similar(integrated.object = myeloid.integrated, assay = 'integrated',layout = 'circular.dendrogram',
#' var.genes.method = 'var', Annotation = 'Annotation', k_clusters = 4, plot = T, return.data = F)
#' ###The distance defined as (1-Pearson correlation coefficient)/2.
#' Find_similar(integrated.object = myeloid.integrated, assay = 'integrated',layout = 'horizontal.dendrogram',var.genes.method = 'var',clustering_distance = 'pearson', Annotation = 'Annotation', plot = T, return.data = F)
#' 
#' @export
Find_similar <- function (integrated.object = NULL, Annotation = "Annotation", 
    var.genes.method = c("var", "pca"), pca.dims = NULL, pca.top.var = 50, 
    assay = "integrated", clustering.method = "complete", clustering_distance = c("cormat", 
        "pearson"), clustering_distance_cols = c("correlation", 
        "euclidean", "binary"), clustering_distance_rows = c("correlation", 
        "euclidean", "binary"), reverse.cor = T, k_clusters = 4, plot = T, return.data = F, 
    layout = c("circular.dendrogram", "horizontal.dendrogram", 
        "heatmap_cc"), circle.text.size = 2.2, circle.line.size = 1, 
    circle.margin = 3, label.color = NULL, dend.color = NULL, 
    heatmap.display.numbers = T) 
{
    var.genes.method <- match.arg(var.genes.method)
    if (var.genes.method == "var") {
        var.genes <- integrated.object@assays[[assay]]@var.features
        if (is.null(var.genes)) stop("var.features doesn't exist.")
    }
    if (var.genes.method == "pca") {
        pca_gene <- integrated.object@reductions$pca@feature.loadings
        if (is.null(pca.dims)) {
            pca.dims <- integrated.object@commands$FindNeighbors.integrated.pca$dims
        }
        pca_gene <- pca_gene[, pca.dims]
        pca_gene <- abs(pca_gene)
        genemat <- apply(pca_gene, 2, function(x) rownames(pca_gene)[order(x, 
            decreasing = T)][1:pca.top.var])
        var.genes <- as.character(genemat)
        if (is.null(var.genes)) stop("PCA doesn't exist.")
    }
    data_vargene <- as.data.frame(Matrix::t(integrated.object@assays[[assay]]@data[var.genes, 
        ]))
    data_vargene$groups <- as.vector(integrated.object@meta.data[, 
        Annotation])
    matm <- reshape2::melt(data_vargene, id.vars = c("groups"), 
        measure.vars = colnames(data_vargene)[-ncol(data_vargene)], 
        variable.name = "var_genes", value.name = "expression")
    colnames(matm) <- c("groups", "var_genes", "expression")
    ave <- tapply(matm$expression, list(matm$groups, matm$var_genes), 
        mean)
    ave <- data.frame(ave, stringsAsFactors = F)
    # remove genes with zero standard deviation.
    sd0_gene <- apply(ave, 2, function(x) length(unique(x)) == 1)
    if (any(sd0_gene)) {
      message(paste0('Warning: remove genes that have same average expression among clusters (the standard deviation is zero): ', colnames(ave)[sd0_gene]))
      ave <- ave[,-which(sd0_gene)]
    }
    matdata <- as.matrix(ave)
    clustering_distance <- match.arg(arg = NULL, clustering_distance)
    clustering_distance_cols <- match.arg(arg = NULL, clustering_distance_cols)
    clustering_distance_rows <- match.arg(arg = NULL, clustering_distance_rows)
    if (clustering_distance == "cormat") {
        cormat <- matrix(0, nrow = nrow(ave), ncol = nrow(ave))
        for (i in 1:nrow(ave)) {
            for (j in i:nrow(ave)) {
                cormat[i, j] <- cor(matdata[i, ], matdata[j, 
                  ])
                cormat[j, i] <- cor(matdata[i, ], matdata[j, 
                  ])
            }
        }
        colnames(cormat) <- rownames(ave)
        rownames(cormat) <- rownames(ave)
        cormat_cc <- cormat
        cormat <- (1 - cormat)/2
        out <- pheatmap::pheatmap(cormat, silent = T, clustering_method = clustering.method, 
            clustering_distance_cols = clustering_distance_cols, 
            clustering_distance_rows = clustering_distance_rows)
    }
    else if (clustering_distance == "pearson") {
        cols.cor <- cor(matdata, use = "pairwise.complete.obs", 
            method = "pearson")
        rows.cor <- cor(t(matdata), use = "pairwise.complete.obs", 
            method = "pearson")
        # reverse.cor: applied method of hierarchical clustering from https://doi.org/10.1016/j.cell.2021.01.010.
        if (reverse.cor) {
          dist.cols <- as.dist( (1 - cols.cor)/2 )
          dist.rows <- as.dist( (1 - rows.cor)/2 )
          # Equal to as.dist(1 - cor)/2
          out <- pheatmap::pheatmap(matdata, silent = T, clustering_method = clustering.method, 
                                    clustering_distance_cols = dist.cols, 
                                    clustering_distance_rows = dist.rows)
          # for heatmap
          cormat_cc <- (1 - rows.cor)/2
        }else{
          out <- pheatmap::pheatmap(rows.cor, silent = T, clustering_method = clustering.method, 
                                    clustering_distance_cols = clustering_distance_cols, 
                                    clustering_distance_rows = clustering_distance_rows)
          # for heatmap
          cormat_cc <- rows.cor
        }
        cormat <- matdata
        
    }
    hc <- out$tree_row
    dend <- as.dendrogram(hc)
    clusters <- cutree(hc, k_clusters)
    if (is.null(label.color)) {
        label_color <- scPalette1(k_clusters)
    }
    else {
        label_color <- label.color
    }
    if (is.null(dend.color)) {
        dend_color <- scPalette1(k_clusters)[unique(clusters[hc$order])]
    }
    else {
        dend_color <- dend.color
    }
    per_label_color <- label_color[clusters[hc$order]]
    dend <- dend %>% dendextend::color_branches(k = k_clusters, 
        col = dend_color)
    ggd1 <- dendextend::as.ggdend(dend)
    ggd1$segments$lwd <- circle.line.size
    labels <- hc$labels[hc$order]
    number_of_type <- length(labels)
    label_data <- data.frame(id = seq(1, number_of_type), individual = labels, 
        per_label_color = per_label_color)
    angle <- 90 - 360 * (label_data$id)/number_of_type
    label_data$hjust <- ifelse(angle < -90, 1.05, -0.05)
    label_data$angle <- ifelse(angle < -90, angle + 180, angle)
    p <- ggplot(ggd1, labels = F) + scale_y_reverse(expand = c(0, 
        0)) + theme_minimal() + theme(axis.text = element_blank(), 
        axis.title = element_blank(), panel.grid = element_blank(), 
        plot.margin = unit(rep(circle.margin, 4), "cm")) + coord_polar(start = 0, 
        clip = "off") + xlim(-0.1, NA) + geom_text(data = label_data, 
        aes(x = id, y = 0, label = individual, hjust = hjust), 
        color = per_label_color, alpha = 1, size = circle.text.size, 
        angle = label_data$angle, inherit.aes = FALSE)
    datalist <- list(cor_mat = cormat, cor_mat_cc = cormat_cc,
                     avedata = matdata, hclust = hc, p = p)
    if (plot) {
        layout <- match.arg(layout)
        if (layout == "heatmap_cc" & clustering_distance == "cormat") {
            pheatmap::pheatmap(datalist[["cor_mat"]], display_numbers = heatmap.display.numbers, 
                clustering_method = clustering.method, clustering_distance_cols = clustering_distance_cols, 
                clustering_distance_rows = clustering_distance_rows)
        }
        if (layout == "heatmap_cc" & clustering_distance == "pearson") {
            pheatmap::pheatmap(datalist[["cor_mat_cc"]], clustering_method = clustering.method, 
                               clustering_distance_cols = clustering_distance_cols, 
                               clustering_distance_rows = clustering_distance_rows)
        }
        if (layout == "horizontal.dendrogram") {
          graphics::plot(datalist[["hclust"]], xlab = "")
        }
        if (layout == "circular.dendrogram") {
            print(datalist[["p"]])
        }
    }
    if (return.data) {
        return(datalist)
    }
}

#' Plot result from Find_similar function.
#' 
#' @examples 
#' myeloid.integrated <- readRDS('/BGFS1/home/huangzh/workspace/project/ovarian_cancer/result/RNA/614.Zhang_Colon_10X_macrophage_12625cells/myeloid.integrated.rds')
#' datalist <- Find_similar(integrated.object = myeloid.integrated, assay = 'integrated',
#' var.genes.method = 'var', Annotation = 'Annotation',k_clusters = 2, plot = F, return.data = T)
#' #correlation heatmap
#' Plot_similar(datalist, layout = 'heatmap',heatmap.display.numbers = T)
#' Plot_similar(datalist, layout = 'heatmap',heatmap.display.numbers = F)
#' #dendrogram
#' Plot_similar(datalist, layout = 'h.dend')
#' # cirlize denrogram
#' Plot_similar(datalist, layout = 'c.dend',k_clusters = 6)
#' 
#' @export
Plot_similar <- function (datalist, k_clusters = 4, layout = c("c.dend", "h.dend", 
    "heatmap","heatmap_cc"), clustering_distance = c("cormat", "pearson"), 
    clustering.method = "complete", clustering_distance_cols = "correlation", 
    clustering_distance_rows = "correlation", circle.text.size = 2.2, 
    circle.line.size = 1, circle.margin = 3, label.color = NULL, 
    per_label_color = NULL, dend.color = NULL, heatmap.display.numbers = T, ...) 
{
    hc <- datalist[["hclust"]]
    dend <- as.dendrogram(hc)
    clusters <- cutree(hc, k_clusters)
    if (is.null(label.color)) {
        label_color <- scPalette1(k_clusters)
    }
    else {
        label_color <- label.color
    }
    if (is.null(dend.color)) {
        dend_color <- scPalette1(k_clusters)[unique(clusters[hc$order])]
    }
    else {
        dend_color <- dend.color[unique(clusters[hc$order])]
    }
    if (is.null(per_label_color)) {
      per_label_color <- label_color[clusters[hc$order]]
    }
    dend <- dend %>% dendextend::color_branches(k = k_clusters, 
        col = dend_color)
    ggd1 <- dendextend::as.ggdend(dend)
    ggd1$segments$lwd <- circle.line.size
    labels <- hc$labels[hc$order]
    number_of_type <- length(labels)
    label_data <- data.frame(id = seq(1, number_of_type), individual = labels, 
        per_label_color = per_label_color)
    angle <- 90 - 360 * (label_data$id)/number_of_type
    label_data$hjust <- ifelse(angle < -90, 1.05, -0.05)
    label_data$angle <- ifelse(angle < -90, angle + 180, angle)
    p <- ggplot(ggd1, labels = F) + scale_y_reverse(expand = c(0, 
        0)) + theme_minimal() + theme(axis.text = element_blank(), 
        axis.title = element_blank(), panel.grid = element_blank(), 
        plot.margin = unit(rep(circle.margin, 4), "cm")) + coord_polar(start = 0, 
        clip = "off") + xlim(-0.1, NA) + geom_text(data = label_data, 
        aes(x = id, y = 0, label = individual, hjust = hjust), 
        color = per_label_color, alpha = 1, size = circle.text.size, 
        angle = label_data$angle, inherit.aes = FALSE)
    datalist <- list(cor_mat = datalist[["cor_mat"]], 
                     cor_mat_cc = datalist[["cor_mat_cc"]],
                     hclust = hc, 
                     p = p)
    layout <- match.arg(layout)
    clustering_distance <- match.arg(arg = NULL, clustering_distance)
    if (layout == "heatmap" & clustering_distance == "cormat") {
        pheatmap::pheatmap(datalist[["cor_mat"]], display_numbers = heatmap.display.numbers, 
            clustering_method = clustering.method, clustering_distance_cols = clustering_distance_cols, 
            clustering_distance_rows = clustering_distance_rows, ...)
    }
    if (layout == "heatmap" & clustering_distance == "pearson") {
        pheatmap::pheatmap(datalist[["avedata"]], clustering_method = clustering.method, 
            clustering_distance_cols = dist.cols, clustering_distance_rows = dist.rows, ...)
    }
    if (layout == "heatmap_cc") {
      pheatmap::pheatmap(datalist[["cor_mat_cc"]], display_numbers = heatmap.display.numbers, 
                         clustering_method = clustering.method, clustering_distance_cols = clustering_distance_cols, 
                         clustering_distance_rows = clustering_distance_rows, ...)
    }
    if (layout == "h.dend") {
      graphics::plot(datalist[["hclust"]], xlab = "")
    }
    if (layout == "c.dend") {
        return(datalist[["p"]])
    }
}