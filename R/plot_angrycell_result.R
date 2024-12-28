#' plot of gene and predicted cell types in GO terms
#'
#' @param object A seurat object.
#' @param active.assay Default assay in seurat object which used to model.
#' @param reduction Select the reduction method.
#' @param angrycell.tab a dataframe of result of angrycell.
#' @param show.all Plot all predicted cell type and marker genes in GO terms.
#' @param return.grob return a list of grob results.
#' @examples
#' anno <- angrycell(pbmc, min.pct = 0.3, show.top = 10, use.definite.markers = F,
#' reduction = 'umap')
#' show_result(pbmc,reduction = 'umap',angrycell.tab = anno,
#' show.all = TRUE,return.grob = FALSE)
#' groblist <- show_result(pbmc,reduction = 'umap',angrycell.tab = anno,
#' show.all = FALSE,return.grob = TRUE)
#' #Because grid.draw only append a new layer on previous plot,
#' #use dev.off() before grid.draw.
#' grid::grid.draw(groblist[['B']])
#' @export
show_result <- function (object = SeuratS4, active.assay = "RNA", reduction = "tsne", 
    angrycell.tab = annotation, show.all = TRUE, return.grob = FALSE, 
    seed = 123) 
{
    if (nrow(angrycell.tab) == 0) {
        stop("Dataframe of Predicted cell types is empty! ")
    }
    if (!any(duplicated(angrycell.tab$Orig_Idents))) {
        stop("Number of candidate predicted cell types in angrycell.tab should be more than 1!\n \n         Parameter show.top in angrycell function should be set to more than 1.")
    }
    set.seed(seed)
    xy <- object@reductions[[reduction]]@cell.embeddings
    gglist <- list()
    for (i in unique(angrycell.tab$Orig_Idents)) {
        plotdata <- angrycell.tab[angrycell.tab$Orig_Idents == 
            i, ]
        if (nrow(plotdata) < 2) {
            warning(paste0("Number of candidate type: ", i, ", must be more than 1."))
            (next)()
        }
        plotdata <- plotdata[!duplicated(plotdata$Celltype_predicted), 
            ]
        genetermlist <- strsplit(as.vector(plotdata$Gene_in_Term), 
            ",")
        names(genetermlist) <- plotdata$Celltype_predicted
        gene_to_plot <- unique(as.vector(unlist(genetermlist)))
        p1 <- highcl(object, i, reduction)
        p1 <- p1 + theme(legend.position = "none") + ggtitle(i)
        if (length(gene_to_plot) > 50) {
            gene_to_plot <- gene_to_plot[1:50]
        }
        exp <- log2(as.data.frame(t(as.matrix(object@assays[[active.assay]]@data[gene_to_plot, 
            ]))) + 1)
        p2 <- fplot1(xy, exp)
        netlist <- list()
        for (n in names(genetermlist)) {
            pos <- gene_to_plot %in% genetermlist[[n]]
            net <- rep(0, length(gene_to_plot))
            net[pos] <- 1
            names(net) <- gene_to_plot
            netlist[[n]] <- net
        }
        netdata <- as.data.frame(do.call(rbind, netlist))
        if (nrow(netdata) == ncol(netdata)) {
            Top_rank <- c(1, rep(0, nrow(netdata) - 1))
            netdata$Top_rank_type <- Top_rank
            gene_to_plot <- c(gene_to_plot, "TOP_RANK")
        }
        bip = network::network(netdata, matrix.type = "bipartite", 
            ignore.eval = FALSE, names.eval = "weights")
        col = c(actor = "gold", event = "green")
        p3 <- GGally::ggnet2(bip, color = "mode", palette = col, 
            label = TRUE, size = c(rep(5, nrow(netdata)), rep(1, 
                length(gene_to_plot)))) + theme(legend.position = "none")
        p <- gridExtra::arrangeGrob(gridExtra::arrangeGrob(p1, 
            p3, ncol = 1, nrow = 2), gridExtra::arrangeGrob(p2, 
            ncol = 1, nrow = 1), heights = c(4, 1), widths = c(1, 
            1))
        gglist[[i]] <- p
        if (show.all) {
            gridExtra::grid.arrange(gridExtra::arrangeGrob(p1, 
                p3, ncol = 1, nrow = 2), gridExtra::arrangeGrob(p2, 
                ncol = 1, nrow = 1), heights = c(4, 1), widths = c(1, 
                1))
        }
    }
    if (return.grob) {
        return(gglist)
    }
}

#' Dotplot for marker gene from the result of angrycell.
#' 
#' @examples 
#' angrydf <- angrycell(pbmc, select.db = c('db_user'), show.top = 2)
#' DotPlot_anno(pbmc, angrydf = angrydf)
#' 
#' @export
DotPlot_anno <- function (object, features = NULL, angrydf = NULL, assay = NULL, 
    auto.order = T, sort.name.per.group.byAnnot = F, human.to.mouse = F, 
    col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
    color = scPalette_heatmap_bar(), group.by = NULL, split.by = NULL, 
    scale.by = "radius", scale.min = NA, scale.max = NA, y.size = 8, 
    x.size = 8, x.vjust = 0.5, x.hjust = 1, x.angle = 90, coord.flip = T, 
    legend.color.name = "Average\nExpression", cols = c("lightgrey","blue"),
    order.method = "DEG") 
{
    if (is.null(features) & !is.null(angrydf)) {
        df_g_to_cell <- lapply(1:nrow(angrydf), function(x) {
            geneset <- as.vector(strsplit(angrydf$Gene_in_Term[x], 
                ",")[[1]])
            dfterm <- data.frame(gene = geneset, cell = angrydf$Celltype_predicted[x], 
                stringsAsFactors = F)
            dfterm
        })
        df_g_to_cell <- unique(do.call(rbind, df_g_to_cell))
        features <- unique(df_g_to_cell$gene)
        if (human.to.mouse) {
            featuredf <- data.frame(V1 = features, stringsAsFactors = F)
            featuredf <- left_join(featuredf, homolog_Human_mouse, 
                by = "V1")
            features <- unique(featuredf$V2)
            features <- features[features %in% rownames(object)]
        }
    }
    if (is.null(x = group.by)) {
        id <- Idents(object = object)
    }
    else {
        id <- object[[group.by, drop = TRUE]]
        Idents(object = object) <- id
    }
    if (auto.order) {
        features <- order_heatmap_gene(object, features, assay = assay, 
            sort.name.per.group.byAnnot = sort.name.per.group.byAnnot, 
            method = order.method,
            gene.Annot = df_g_to_cell[, c("cell", "gene")])
    }
    if (!is.null(assay)) {
        DefaultAssay(object = object) <- assay
    }
    scale.func <- switch(EXPR = scale.by, size = scale_size, 
        radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
    data.features <- FetchData(object = object, vars = features)
    data.features$id <- id
    if (!is.factor(x = data.features$id)) {
        data.features$id <- factor(x = data.features$id)
    }
    id.levels <- levels(x = data.features$id)
    data.features$id <- as.vector(x = data.features$id)
    if (!is.null(x = split.by)) {
        splits <- object[[split.by, drop = TRUE]]
        if (length(x = unique(x = splits)) > length(x = cols)) {
            stop("Not enought colors for the number of groups")
        }
        cols <- cols[1:length(x = unique(x = splits))]
        names(x = cols) <- unique(x = splits)
        data.features$id <- paste(data.features$id, splits, sep = "_")
        unique.splits <- unique(x = splits)
        id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
            "_", rep(x = unique(x = splits), times = length(x = id.levels)))
    }
    data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
        data.use <- data.features[data.features$id == ident, 
            1:(ncol(x = data.features) - 1), drop = FALSE]
        avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
            return(mean(x = expm1(x = x)))
        })
        pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, 
            threshold = 0)
        return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    })
    names(x = data.plot) <- unique(x = data.features$id)
    data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
        data.use <- as.data.frame(x = data.plot[[x]])
        data.use$features.plot <- rownames(x = data.use)
        data.use$id <- x
        return(data.use)
    })
    data.plot <- do.call(what = "rbind", args = data.plot)
    if (!is.null(x = id.levels)) {
        data.plot$id <- factor(x = data.plot$id, levels = id.levels)
    }
    avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
        FUN = function(x) {
            data.use <- data.plot[data.plot$features.plot == 
                x, "avg.exp"]
            if (length(x = data.use) > 1) {
                data.use <- scale(x = data.use)
            }
            data.use <- MinMax(data = data.use, min = col.min, 
                max = col.max)
            return(data.use)
        })
    avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
    if (!is.null(x = split.by)) {
        avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
            breaks = 20))
    }
    data.plot$avg.exp.scaled <- avg.exp.scaled
    data.plot$features.plot <- factor(x = data.plot$features.plot, 
        levels = rev(x = features))
    data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
    data.plot$pct.exp <- data.plot$pct.exp * 100
    if (!is.null(x = split.by)) {
        splits.use <- vapply(X = strsplit(x = as.character(x = data.plot$id), 
            split = "_"), FUN = "[[", FUN.VALUE = character(length = 1L), 
            2)
        data.plot$colors <- mapply(FUN = function(color, value) {
            return(colorRampPalette(colors = c("grey", color))(20)[value])
        }, color = cols[splits.use], value = avg.exp.scaled)
    }
    color.by <- ifelse(test = is.null(x = split.by), yes = "avg.exp.scaled", 
        no = "colors")
    if (!is.na(x = scale.min)) {
        data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
    }
    if (!is.na(x = scale.max)) {
        data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
    }
    if (!is.null(angrydf)) {
        df_g_to_cell <- lapply(1:nrow(angrydf), function(x) {
            geneset <- as.vector(strsplit(angrydf$Gene_in_Term[x], 
                ",")[[1]])
            if (human.to.mouse) {
                featuredf <- data.frame(V1 = geneset, stringsAsFactors = F)
                featuredf <- left_join(featuredf, homolog_Human_mouse, 
                  by = "V1")
                geneset <- unique(featuredf$V2)
                geneset <- geneset[geneset %in% rownames(object)]
            }
            dfterm <- data.frame(gene = geneset, cell = angrydf$Celltype_predicted[x], 
                stringsAsFactors = F)
            dfterm
        })
        df_g_to_cell <- unique(do.call(rbind, df_g_to_cell))
        list_g_to_cell <- sapply(unique(df_g_to_cell$gene), function(x) {
            cellset <- gsub(" ", "", df_g_to_cell$cell[df_g_to_cell$gene == 
                x])
            cells <- paste(cellset, collapse = ", ")
        })
        features_annodf <- data.frame(features.plot = unique(df_g_to_cell$gene), 
            anno = list_g_to_cell, stringsAsFactors = F)
        leveldf <- data.frame(features.plot = rev(features), 
            index = 1:length(features), stringsAsFactors = F)
        data.plot <- left_join(data.plot, leveldf, by = "features.plot")
        data.plot <- left_join(data.plot, features_annodf, by = "features.plot")
        pos <- !is.na(data.plot$anno)
        data.plot$features.plot[pos] <- paste("(", data.plot$anno[pos], 
            ") ", data.plot$features.plot[pos])
        leveldf2 <- unique(data.plot[, c("features.plot", "index")])
        leveldf2 <- leveldf2[order(leveldf2$index), ]
        data.plot$features.plot <- factor(data.plot$features.plot, 
            levels = leveldf2$features.plot)
    }
    plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", 
        y = "id")) + geom_point(mapping = aes_string(size = "pct.exp", 
        color = color.by)) + geom_point(mapping = aes_string(size = "pct.exp"), 
        color = "black", pch = 21) + scale.func(range = c(0, 
        dot.scale), limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) + 
        labs(x = "Features", y = ifelse(test = is.null(x = split.by), 
            yes = "Identity", no = "Split Identity")) + cowplot::theme_cowplot()
    plot <- plot + scale_color_gradientn(name = legend.color.name, 
        colors = rev(color)) + theme_bw() + xlab("") + ylab("") + 
        theme(axis.text.x = element_text(size = x.size, color = "black", 
            vjust = x.vjust, hjust = x.hjust, angle = x.angle), 
            axis.text.y = element_text(size = y.size, color = "black"), 
            panel.grid.major = element_line(size = 0.1, linetype = "dotted", 
                colour = "grey"))
    if (coord.flip) {
        plot <- plot + coord_flip()
    }
    return(plot)
}

#' Show predicted cell types in tree。
#'
#' @export
List_predicted_celltypes <- function(angrycell.tab, outfile = NULL, type = "text",
                                     root = "Cell types"){
  prdcell <- angrycell.tab$Celltype_predicted[!duplicated(angrycell.tab$Orig_Idents)]
  list1 <- mindrlist1[names(mindrlist1) %in% prdcell]
  list2 <- mindrlist2[names(mindrlist2) %in% prdcell]
  list2 <- list()
  for (t in names(mindrlist2)) {
    value <- mindrlist2[[t]][mindrlist2[[t]] %in% prdcell]
    if (length(value) == 0) {
      next()
    }
    list2[[t]] <- value
  }
  list1 <- list()
  for (t in names(mindrlist1)) {
    value <- mindrlist1[[t]][mindrlist1[[t]] %in% names(list2)]
    if (length(value) == 0) {
      next()
    }
    list1[[t]] <- value
  }
  
  input <- c()
  for (i in names(list1)) {
    type1 <- paste0('# ',i)
    input <- c(input, type1)
    if (length(list1[[i]]) > 0) {
      type2list <-list1[[i]]
      for (j in type2list) {
        input <- c(input,  paste0('## ', j))
        if (length(list2[[j]]) > 0) {
          input <- c(input,  paste0('### ', list2[[j]]))
        }
      }
    }
  }
  mindr::mm(from = input, to = outfile, type = type, root = root)
}
#' Plot adjust p value of cell types.
#'
#' @param candidate Genes to test.
#' @param N Number of cells.
#' @param db Marker genes database for usage.
#' @return ggplot object.
#' @examples
#' anno <- angrycell(pbmc)
#' pbmc@@meta.data$Orig_Idents <- Idents(pbmc)
#' pbmc@@meta.data$Cellname <- rownames(pbmc@@meta.data)
#' pbmc@@meta.data <- inner_join(pbmc@@meta.data,anno)
#' rownames(pbmc@@meta.data) <- pbmc@@meta.data$Cellname
#' plot_pvalue(pbmc)
#' @export
plot_pvalue <- function (objs4, reduction = "umap") 
{
  mycol <- RColorBrewer::brewer.pal(11, "RdYlBu")
  objs4@meta.data$Orig_Idents <- as.vector(gsub(" ", "_", objs4@meta.data$Orig_Idents))
  objs4@meta.data$id_p <- paste0(objs4@meta.data$Orig_Idents, 
                                 "_", objs4@meta.data$BH_adj_Pvalue)
  Idents(objs4) <- objs4@meta.data$id_p
  p <- DimPlot(objs4, label = T, repel = T) + NoLegend()
  xy <- objs4@reductions[[reduction]]@cell.embeddings
  xy <- as.data.frame(xy)
  xy$BH_adj_Pvalue <- objs4@meta.data$BH_adj_Pvalue
  p1 <- ggplot(xy, aes(x = xy[, 1], y = xy[, 2], colour = BH_adj_Pvalue)) + 
    geom_point(size = 0.01) + theme_bw() + theme(panel.grid.major = element_line(colour = NA), 
                                                 panel.grid.minor = element_line(colour = NA))
  p1 <- p1 + scale_colour_gradient2(low = mycol[5], mid = mycol[2], 
                                    high = mycol[1], midpoint = 0.5) + xlab(colnames(xy)[1]) + 
    ylab(colnames(xy)[2])
  return(p + p1)
}