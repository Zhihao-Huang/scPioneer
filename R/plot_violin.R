#' scater plot of gene expression.
#'
#' @return A ggplot object.
#' @param data a normalized expression matrix(gene x cell).
#' @param gene Genes to plot.
#' @param ident Cell annotation. Idents(seurat_object)
#' @param facet.toward 'col' or 'row'.
#' @param strip.text.size Gene name text size.
#' @examples
#' data('pbmc')
#' marker_to_plot <- c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A")
#' exp <- as.data.frame(t(as.matrix(pbmc@@assays$RNA@@data[marker_to_plot,, drop=FALSE])))
#' vioplot(exp, marker_to_plot, Idents(pbmc))
#' vioplot(exp, marker_to_plot, Idents(pbmc),fill.expression = T,add.ave.point = T,add.line = T,
#' jitter = T,add.box = T,legend.name = 'Average expression')
#' @export
vioplot <- function (data, paths, ident, select.celltypes = NULL, facet.toward = "row", 
    do.log = F, strip.text.size = 8, strip.text.angle = NULL, 
    strip.text.hjust = NULL, strip.text.vjust = NULL, panel.spacing = 0.5, 
    clip = "on", color.use = NULL, jitter = F, cell.order = NULL, 
    add.ave.point = T, ave.point.size = 1.5, add.line = T, line.size = 0.5, 
    add.box = F, strip.text.switch = NULL, border.color = "black", 
    border.size = 0.5, fill.expression = F, fill.color = c("#FEE090", 
        "#D73027", "#A50026"), mid.color.point = "auto", legend.position = "top", 
    legend.name = "", jitter.size = 0.1, text.x.size = 8, text.y.size = 8, 
    hjust = 0, vjust = 0, legend.scale = NULL, legend.title.size = 8, 
    legend.text.size = 8, legend.key.size = 1, axis.text.x.angle = NULL, 
    remove.axis.x = F, remove.axis.y = F) 
{
    sig_path <- as.data.frame(data[, paths, drop = FALSE])
    if (do.log) {
        sig_path <- log2(sig_path + 1)
    }
    sig_path$groups <- ident
    if (is.null(levels(sig_path$groups))) {
        sig_path$groups <- factor(sig_path$groups, levels = sort(unique(sig_path$groups)))
    }
    is_all_numeric <- sum(is.na(suppressWarnings(as.numeric(levels(sig_path$groups))))) == 
        0
    if (is.null(cell.order) & is_all_numeric) {
        message("Warning: add 'c_' in front of numeric cluster names to keep cell order.")
        newclusters <- paste0("c_", as.vector(sig_path$groups))
        newlevels <- paste0("c_", sort(as.numeric(levels(sig_path$groups))))
        ident <- factor(newclusters, levels = newlevels)
        sig_path$groups <- ident
    }
    if (!is.null(select.celltypes)) {
        pos <- sig_path$groups %in% select.celltypes
        if (sum(pos) == 0) {
            stop("All types specified by select.celtypes are not in ident.")
        }
        pos2 <- select.celltypes %in% sig_path$groups
        if (sum(pos2) != length(pos2)) {
            message(paste0("Warning: cell types that were not in ident: ", 
                paste(select.celltypes[!pos2], collapse = ", ")))
        }
        sig_path <- sig_path[pos, ]
        sig_path$groups <- factor(sig_path$groups, levels = select.celltypes[pos2])
    }
    matm <- melt(sig_path, id.vars = c("groups"), measure.vars = paths, 
        variable.name = "pathway", value.name = "TPM")
    colnames(matm) <- c("groups", "pathway", "TPM")
    ave <- tapply(matm$TPM, list(matm$groups, matm$pathway), 
        mean)
    ave <- as.data.frame(ave)
    ave$groups <- rownames(ave)
    ave <- reshape2::melt(ave, id.vars = c("groups"), measure.vars = paths, 
        variable.name = "pathway", value.name = "ave", factorsAsStrings = T)
    matm <- suppressMessages(inner_join(matm, ave))
    if (!is.null(cell.order)) {
        matm$groups <- factor(matm$groups, levels = cell.order)
    }
    color_palette = c("#00A087FF", "#4DBBD5FF", "#E64B35FF", 
        "#3C5488FF", "#F38400", "#A1CAF1", "#BE0032", "#C2B280", 
        "#848482", "#008856", "#E68FAC", "#0067A5", "#604E97", 
        "#F6A600", "#B3446C", "#DCD300", "#882D17", "#8DB600", 
        "#654522", "#E25822", "#2B3D26", "#848482", "#008856", 
        "#E68FAC", "#0067A5", "#604E97", "#F6A600", "#B3446C", 
        "#DCD300", "#882D17")
    if (fill.expression) {
        if (mid.color.point == "auto") {
            mid.color <- (max(matm$ave) - min(matm$ave))/2
        }
        else {
            mid.color <- mid.color.point
        }
        p <- ggplot(matm, aes(x = groups, y = TPM, fill = ave))
        p <- p + geom_violin(trim = T, scale = "width", color = border.color, 
            size = border.size) + theme_bw()
        if (is.null(legend.scale)) {
            p <- p + scale_fill_gradient2(name = legend.name, 
                low = fill.color[1], mid = fill.color[2], high = fill.color[3], 
                midpoint = mid.color)
        }
        else if ("auto" %in% legend.scale) {
            legend.scale <- c(0, max(matm$ave))
            p <- p + scale_fill_gradientn(colours = fill.color, 
                values = scales::rescale(legend.scale), guide = "colorbar", 
                limits = legend.scale, name = legend.name)
        }
        else {
            p <- p + scale_fill_gradientn(colours = fill.color, 
                values = scales::rescale(legend.scale), guide = "colorbar", 
                limits = legend.scale, name = legend.name)
        }
    }
    else {
        p <- ggplot(matm, aes(x = groups, y = TPM, fill = groups))
        p <- p + geom_violin(scale = "width", color = border.color, 
            size = border.size) + theme_bw()
        if (is.null(color.use) & length(unique(matm$groups)) < 
            31) {
            p <- p + scale_fill_manual(values = color_palette)
        }
        else {
            p <- p + scale_fill_manual(values = color.use)
        }
        legend.position <- "none"
    }
    if (jitter) {
        p <- p + geom_jitter(size = jitter.size)
    }
    if (add.ave.point) {
        p <- p + geom_point(data = matm, aes(x = groups, y = ave, 
            group = 1), color = "Black", shape = 15, size = ave.point.size)
    }
    if (add.line) {
        p <- p + geom_line(data = matm, aes(x = groups, y = ave, 
            group = 1), color = "Black", size = line.size, linetype = "dashed")
    }
    if (add.box) {
        p <- p + geom_boxplot(width = 0.1, fill = "white", color = "Black", 
            outlier.shape = NA)
    }
    if (facet.toward == "row") {
        axis.text.x.angle <- ifelse(is.null(axis.text.x.angle), 
            90, axis.text.x.angle)
        p <- p + facet_grid(rows = vars(pathway), scales = "free", 
            switch = strip.text.switch) + theme(legend.position = legend.position, 
            panel.grid.major = element_line(colour = NA), panel.grid.minor = element_line(colour = NA)) + 
            theme(axis.text.x = element_text(angle = axis.text.x.angle, 
                size = text.x.size, hjust = hjust, vjust = vjust), 
                axis.text.y = element_text(size = text.y.size))
    }
    else if (facet.toward == "col") {
        axis.text.x.angle <- ifelse(is.null(axis.text.x.angle), 
            270, axis.text.x.angle)
        p <- p + facet_grid(cols = vars(pathway), scales = "free", 
            switch = strip.text.switch) + theme(legend.position = legend.position, 
            panel.grid.major = element_line(colour = NA), panel.grid.minor = element_line(colour = NA)) + 
            theme(axis.text.x = element_text(angle = axis.text.x.angle, 
                size = text.x.size, hjust = hjust, vjust = vjust), 
                axis.text.y = element_text(size = text.y.size))
        p <- p + coord_flip()
    }
    if (!is.null(strip.text.switch)) {
        p <- p + scale_y_continuous(position = "right")
    }
    if (remove.axis.x) {
        p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
            axis.ticks.x = element_blank())
    }
    if (remove.axis.y) {
        p <- p + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), 
            axis.ticks.y = element_blank())
    }
    p <- p + xlab("") + ylab("")
    p <- p + theme(legend.title = element_text(size = legend.title.size), 
        legend.text = element_text(size = legend.text.size), 
        legend.key.size = unit(legend.key.size, "lines"), strip.text = element_text(size = strip.text.size, 
            angle = strip.text.angle, hjust = strip.text.hjust, 
            vjust = strip.text.vjust), strip.background = element_blank(), 
        strip.placement = "outside", panel.spacing = unit(panel.spacing, 
            "lines"))
    p <- p + guides(color = FALSE)
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
#' violin plot of gene expression.
#'
#' See ?vioplot for R documentation.
#' 
#' @param object A Seurat object.
#' @param paths Genes to plot.
#' @param ident Cell types or other meta data.
#' @param ... Arguments passed to vioplot.
#' 
#' @return A ggplot object.
#' 
#' @examples
#' marker_to_plot <- c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A")
#' vioplot2(pbmc, marker_to_plot)
#' vioplot2(pbmc, marker_to_plot, fill.expression = F,add.ave.point = F,add.line = F,jitter = F,add.box = F)
#' vioplot2(pbmc, marker_to_plot, fill.expression = F,add.ave.point = F,add.line = F,jitter = F,add.box = T)
#' vioplot2(pbmc, marker_to_plot, fill.expression = F,add.ave.point = T,add.line = T,jitter = F,add.box = T)
#' vioplot2(pbmc, marker_to_plot, fill.expression = F,add.ave.point = T,add.line = T,jitter = T,add.box = T)
#' vioplot2(pbmc, marker_to_plot, fill.expression = T,add.ave.point = T,add.line = T,jitter = F,add.box = F)
#' vioplot2(pbmc, marker_to_plot, fill.expression = F,add.ave.point = T,add.line = T,jitter = T,add.box = T,
#' facet.toward = 'col')
#' 
#' @export
vioplot2 <- function (object, paths, ident = Idents(object), assay = "RNA", 
          select.celltypes = NULL, facet.toward = "row", do.log = F, 
          strip.text.size = 8, strip.text.angle = NULL, strip.text.hjust = NULL, 
          strip.text.vjust = NULL, strip.text.switch = NULL, border.color = "black", 
          border.size = 0.5, color.use = NULL, jitter = F, cell.order = NULL, 
          add.ave.point = F, add.line = F, line.size = 0.5, add.box = F, 
          fill.expression = F, fill.color = c("#FEE090", "#D73027", "#A50026"), 
          mid.color.point = "auto",
          legend.position = "top", 
          legend.name = "", jitter.size = 0.1, text.x.size = 8, text.y.size = 8, 
          hjust = 1, vjust = 1, legend.scale = NULL, legend.title.size = 8, 
          axis.text.x.angle = 45, panel.spacing = 0.5, clip = "on", 
          remove.axis.x = F, remove.axis.y = F, ...) 
{
  data <- as.data.frame(t(as.matrix(GetAssayData(object, assay = assay)[paths, 
                                                                , drop = FALSE])))
  p <- vioplot(data, paths, ident = ident, select.celltypes = select.celltypes, 
               facet.toward = facet.toward, do.log = do.log, strip.text.size = strip.text.size, 
               strip.text.angle = strip.text.angle, strip.text.hjust = strip.text.hjust, 
               strip.text.vjust = strip.text.vjust, strip.text.switch = strip.text.switch, 
               border.color = border.color, border.size = border.size, 
               color.use = color.use, jitter = jitter, cell.order = cell.order, 
               add.ave.point = add.ave.point, add.line = add.line, line.size = line.size, 
               add.box = add.box, fill.expression = fill.expression, 
               fill.color = fill.color, mid.color.point = mid.color.point, 
               legend.position = legend.position, legend.name = legend.name, 
               jitter.size = jitter.size, text.x.size = text.x.size, 
               text.y.size = text.y.size, hjust = hjust, vjust = vjust, 
               legend.scale = legend.scale, legend.title.size = legend.title.size, 
               axis.text.x.angle = axis.text.x.angle, panel.spacing = panel.spacing, 
               clip = clip, remove.axis.x = remove.axis.x, remove.axis.y = remove.axis.y, 
               ...)
  return(p)
}
#' Violin plot of average expression of gene set.
#' @return A ggplot object.
#' @param pathlist A list includes several gene sets.
#' @export
vioplot3 <- function (object, pathlist, ident = Idents(object), assay = "RNA", 
          ...) 
{
  datalist <- list()
  names(pathlist) <- gsub(" ", ".", names(pathlist))
  for (i in names(pathlist)) {
    pos <- pathlist[[i]] %in% rownames(object)
    message(paste0("Genes of ", i, " not in matrix: ", pathlist[[i]][!pos]))
    genes <- pathlist[[i]][pos]
    data <- as.data.frame(t(as.matrix(GetAssayData(object, assay = assay)[genes, 
                                                                  , drop = FALSE])))
    datalist[[i]] <- apply(data, 1, mean)
  }
  data <- data.frame(datalist, stringsAsFactors = F)
  p <- vioplot(data, names(pathlist), ident = ident, ...)
  return(p)
}

#' violin plot for one celltype
#' 
#' @export
vioplot_one_celltype <- function (object, paths, ident = Idents(object), assay = "RNA", 
    facet.toward = "row", do.log = F, strip.text.size = 8, strip.text.angle = NULL, 
    strip.text.hjust = NULL, strip.text.vjust = NULL, strip.text.switch = NULL, 
    border.color = "black", border.size = 0.5, color.use = NULL, 
    jitter = F, cell.order = NULL, add.ave.point = T, add.line = T, 
    line.size = 0.5, add.box = F, fill.expression = F, fill.color = c("#FEE090", 
        "#D73027", "#A50026"), mid.color.point = "auto", legend.position = "top", 
    legend.name = "", jitter.size = 0.1, text.x.size = 8, text.y.size = 8, 
    hjust = 0, vjust = 0, legend.scale = NULL, legend.title.size = 8, 
    axis.text.x.angle = NULL, panel.spacing = 0.5, clip = "on", 
    remove.axis.x = F, remove.axis.y = F, ave.point.size = 1.5, 
    legend.key.size = 1, legend.text.size = 8) 
{
    data <- as.data.frame(t(as.matrix(GetAssayData(object, assay = assay)[paths, 
        , drop = FALSE])))
    sig_path <- as.data.frame(data[, paths, drop = FALSE])
    if (do.log) {
        sig_path <- log2(sig_path + 1)
    }
    sig_path$groups <- ident
    matm <- melt(sig_path, id.vars = c("groups"), measure.vars = paths, 
        variable.name = "pathway", value.name = "TPM")
    colnames(matm) <- c("groups", "pathway", "TPM")
    ave <- tapply(matm$TPM, list(matm$groups, matm$pathway), 
        mean)
    ave <- as.data.frame(ave)
    ave$groups <- rownames(ave)
    ave <- reshape2::melt(ave, id.vars = c("groups"), measure.vars = paths, 
        variable.name = "pathway", value.name = "ave", factorsAsStrings = T)
    matm <- suppressMessages(inner_join(matm, ave))
    if (!is.null(cell.order)) {
        matm$groups <- factor(matm$groups, levels = cell.order)
    }
    color_palette = c("#00A087FF", "#4DBBD5FF", "#E64B35FF", 
        "#3C5488FF", "#F38400", "#A1CAF1", "#BE0032", "#C2B280", 
        "#848482", "#008856", "#E68FAC", "#0067A5", "#604E97", 
        "#F6A600", "#B3446C", "#DCD300", "#882D17", "#8DB600", 
        "#654522", "#E25822", "#2B3D26", "#848482", "#008856", 
        "#E68FAC", "#0067A5", "#604E97", "#F6A600", "#B3446C", 
        "#DCD300", "#882D17")
    if (fill.expression) {
        if (mid.color.point == "auto") {
            mid.color <- (max(matm$ave) - min(matm$ave))/2
        }
        else {
            mid.color <- mid.color.point
        }
        p <- ggplot(matm, aes(x = pathway, y = TPM, fill = ave))
        p <- p + geom_violin(trim = T, scale = "width", color = border.color, 
            size = border.size) + theme_bw()
        if (is.null(legend.scale)) {
            p <- p + scale_fill_gradient2(name = legend.name, 
                low = fill.color[1], mid = fill.color[2], high = fill.color[3], 
                midpoint = mid.color)
        }
        else if ("auto" %in% legend.scale) {
            legend.scale <- c(0, max(matm$ave))
            p <- p + scale_fill_gradientn(colours = fill.color, 
                values = scales::rescale(legend.scale), guide = "colorbar", 
                limits = legend.scale, name = legend.name)
        }
        else {
            p <- p + scale_fill_gradientn(colours = fill.color, 
                values = scales::rescale(legend.scale), guide = "colorbar", 
                limits = legend.scale, name = legend.name)
        }
    }
    else {
        p <- ggplot(matm, aes(x = pathway, y = TPM, fill = pathway))
        p <- p + geom_violin(scale = "width", color = border.color, 
            size = border.size) + theme_bw()
        if (is.null(color.use) & length(unique(matm$pathway)) < 
            31) {
            p <- p + scale_fill_manual(values = color_palette)
        }
        else {
            p <- p + scale_fill_manual(values = color.use)
        }
        legend.position <- "none"
    }
    if (jitter) {
        p <- p + geom_jitter(size = jitter.size)
    }
    if (add.ave.point) {
        p <- p + geom_point(data = matm, aes(x = pathway, y = ave, 
            group = 1), color = "Black", shape = 15, size = ave.point.size)
    }
    if (add.line) {
        p <- p + geom_line(data = matm, aes(x = pathway, y = ave, 
            group = 1), color = "Black", size = line.size, linetype = "dashed")
    }
    if (add.box) {
        p <- p + geom_boxplot(width = 0.1, fill = "white", color = "Black", 
            outlier.shape = NA)
    }
    if (facet.toward == "row") {
        axis.text.x.angle <- ifelse(is.null(axis.text.x.angle), 
            90, axis.text.x.angle)
        p <- p + facet_grid(rows = vars(groups), scales = "free", 
            switch = strip.text.switch) + theme(legend.position = legend.position, 
            panel.grid.major = element_line(colour = NA),
            panel.grid.minor = element_line(colour = NA)) + 
            theme(axis.text.x = element_text(angle = axis.text.x.angle, 
                size = text.x.size, hjust = hjust, vjust = vjust), 
                axis.text.y = element_text(size = text.y.size))
    }
    else if (facet.toward == "col") {
        axis.text.x.angle <- ifelse(is.null(axis.text.x.angle), 
            270, axis.text.x.angle)
        p <- p + facet_grid(cols = vars(groups), scales = "free", 
            switch = strip.text.switch) + theme(legend.position = legend.position, 
            panel.grid.major = element_line(colour = NA), 
            panel.grid.minor = element_line(colour = NA)) + 
            theme(axis.text.x = element_text(angle = axis.text.x.angle, 
                size = text.x.size, hjust = hjust, vjust = vjust), 
                axis.text.y = element_text(size = text.y.size))
        p <- p + coord_flip()
    }
    if (!is.null(strip.text.switch)) {
        p <- p + scale_y_continuous(position = "right")
    }
    if (remove.axis.x) {
        p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
            axis.ticks.x = element_blank())
    }
    if (remove.axis.y) {
        p <- p + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), 
            axis.ticks.y = element_blank())
    }
    p <- p + xlab("") + ylab("")
    p <- p + theme(legend.title = element_text(size = legend.title.size), 
        legend.text = element_text(size = legend.text.size), 
        legend.key.size = unit(legend.key.size, "lines"), 
        strip.text = element_text(size = strip.text.size, 
            angle = strip.text.angle, hjust = strip.text.hjust, 
            vjust = strip.text.vjust),
        strip.background = element_blank(), 
        strip.placement = "outside", panel.spacing = unit(panel.spacing, 
            "lines"))
    p <- p + guides(color = FALSE)
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
#' Violin plot of one gene in different groups for each cell type.
#' 
#' @param celltype A vector contains cell types (same as Idents(object));
#' @param group A vector contains group names (group of cell types).
#' @param gene Character. One gene name.
#' 
#' @export
violin_facet_one_gene <- function (object, celltype, group, gene, assay = "RNA", facet.toward = "row", 
    strip.text.size = 8, strip.text.angle = NULL, strip.text.hjust = NULL, 
    strip.text.vjust = NULL, panel.spacing = 0.5, clip = "on", 
    border.color = "black", border.size = 0.5, color.use = NULL, 
    jitter = F, celltype.order = NULL, group.order = NULL, add.ave.point = T, 
    ave.point.size = 1.5, add.line = T, line.size = 0.5, add.box = F, 
    strip.text.switch = NULL, fill.expression = F, fill.color = brewer.pal(n = 9, 
        name = "YlOrRd"), mid.color.point = "auto", legend.position = "top", 
    legend.name = "", jitter.size = 0.1, text.x.size = 8, text.y.size = 8, 
    hjust = 0, vjust = 0, legend.scale = NULL, legend.title.size = 8, 
    axis.text.x.angle = NULL, remove.axis.x = F, remove.axis.y = F) 
{
    data <- as.data.frame(t(as.matrix(GetAssayData(object, assay = assay)[gene, 
        , drop = FALSE])))
    data$celltype <- celltype
    data$group <- group
    ave <- tapply(data[,gene], list(data$group, data$celltype), 
        mean)
    ave <- as.data.frame(ave)
    ave$group <- rownames(ave)
    ave <- reshape2::melt(ave, id.vars = c("group"), measure.vars = unique(celltype), 
        variable.name = "celltype", value.name = "ave", factorsAsStrings = T)
    data <- inner_join(data, ave, by = c("celltype", "group"))
    colnames(data)[1] <- "gene"
    if (!is.null(group.order)) {
        data$group <- factor(data$group, levels = group.order)
    }
    if (!is.null(celltype.order)) {
        data$celltype <- factor(data$celltype, levels = celltype.order)
    }
    color_palette = c("#00A087FF", "#4DBBD5FF", "#E64B35FF", 
        "#3C5488FF", "#F38400", "#A1CAF1", "#BE0032", "#C2B280", 
        "#848482", "#008856", "#E68FAC", "#0067A5", "#604E97", 
        "#F6A600", "#B3446C", "#DCD300", "#882D17", "#8DB600", 
        "#654522", "#E25822", "#2B3D26", "#848482", "#008856", 
        "#E68FAC", "#0067A5", "#604E97", "#F6A600", "#B3446C", 
        "#DCD300", "#882D17")
    if (fill.expression) {
        if (mid.color.point == "auto") {
            mid.color <- (max(data$ave) - min(data$ave))/2
        }
        else {
            mid.color <- mid.color.point
        }
        p <- ggplot(data, aes(x = group, y = gene, fill = ave))
        p <- p + geom_violin(trim = T, scale = "width", color = border.color, 
            size = border.size) + theme_bw()
        if (is.null(legend.scale)) {
            p <- p + scale_fill_gradient2(name = legend.name, 
                low = fill.color[1], mid = fill.color[2], high = fill.color[3], 
                midpoint = mid.color)
        }
        else if (legend.scale == "auto") {
            legend.scale <- c(0, max(data$ave))
            p <- p + scale_fill_gradientn(colours = fill.color, 
                values = scales::rescale(legend.scale), guide = "colorbar", 
                limits = legend.scale, name = legend.name)
        }
        else {
            p <- p + scale_fill_gradientn(colours = fill.color, 
                values = scales::rescale(legend.scale), guide = "colorbar", 
                limits = legend.scale, name = legend.name)
        }
    }
    else {
        p <- ggplot(data, aes(x = group, y = gene, color = group))
        p <- p + geom_violin(fill = NA, scale = "width", color = border.color, 
            size = border.size) + theme_bw()
        if (is.null(color.use) & length(unique(data$group)) < 
            31) {
            p <- p + scale_color_manual(values = color_palette)
        }
        else {
            p <- p + scale_color_manual(values = color.use)
        }
        legend.position <- "none"
    }
    if (jitter) {
        p <- p + geom_jitter(size = jitter.size)
    }
    if (add.ave.point) {
        p <- p + geom_point(data = data, aes(x = group, y = ave, 
            group = 1), color = "Black", shape = 15, size = ave.point.size)
    }
    if (add.line) {
        p <- p + geom_line(data = data, aes(x = group, y = ave, 
            group = 1), color = "Black", size = line.size, linetype = "dashed")
    }
    if (add.box) {
        p <- p + geom_boxplot(width = 0.1, fill = "white", color = "Black", 
            outlier.shape = NA)
    }
    if (facet.toward == "row") {
        axis.text.x.angle <- ifelse(is.null(axis.text.x.angle), 
            90, axis.text.x.angle)
        p <- p + facet_grid(rows = vars(celltype), scales = "free", 
            switch = strip.text.switch) + theme(legend.position = legend.position, 
            panel.grid.major = element_line(colour = NA), panel.grid.minor = element_line(colour = NA)) + 
            theme(axis.text.x = element_text(angle = axis.text.x.angle, 
                size = text.x.size, hjust = hjust, vjust = vjust), 
                axis.text.y = element_text(size = text.y.size))
        p <- p + coord_cartesian(clip = clip)
    }
    else if (facet.toward == "col") {
        axis.text.x.angle <- ifelse(is.null(axis.text.x.angle), 
            270, axis.text.x.angle)
        p <- p + coord_flip(clip = clip)
        p <- p + facet_grid(cols = vars(celltype), scales = "free", 
            switch = strip.text.switch) + theme(legend.position = legend.position, 
            panel.grid.major = element_line(colour = NA), panel.grid.minor = element_line(colour = NA)) + 
            theme(axis.text.x = element_text(angle = axis.text.x.angle, 
                size = text.x.size, hjust = hjust, vjust = vjust), 
                axis.text.y = element_text(size = text.y.size))
    }
    if (!is.null(strip.text.switch)) {
        p <- p + scale_y_continuous(position = "right")
    }
    if (remove.axis.x) {
        p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
            axis.ticks.x = element_blank())
    }
    if (remove.axis.y) {
        p <- p + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), 
            axis.ticks.y = element_blank())
    }
    p <- p + xlab("") + ylab("")
    p <- p + theme(legend.title = element_text(size = legend.title.size), 
        strip.text = element_text(size = strip.text.size, angle = strip.text.angle, 
            hjust = strip.text.hjust, vjust = strip.text.vjust), 
        strip.background = element_blank(), strip.placement = "outside", 
        panel.spacing = unit(panel.spacing, "lines"))
    p <- p + guides(color = FALSE)
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