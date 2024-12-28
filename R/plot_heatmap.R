#' heatmap of average expression
#' @export
aveheatmap <- function(object, paths, ident = Idents(object),  assay = 'RNA',
                       anno_col = NULL, anno_row = NULL) {
  data <- as.data.frame(t(as.matrix(object@assays[[assay]]@data[paths,, drop=FALSE])))
  sig_path <- as.data.frame(data[,paths, drop=FALSE])
  sig_path <- log2(sig_path + 1)
  sig_path$groups <- ident
  matm <- melt(sig_path, id.vars = c("groups"), measure.vars = paths,
               variable.name = "pathway", value.name = "TPM")
  # attach(matm)
  colnames(matm) <- c("groups", "pathway", "TPM")
  ave <- tapply(matm$TPM, list(matm$groups, matm$pathway), mean)
  ave <- as.data.frame(ave)
  ave <- ave[levels(ident),]
  color_palette = c("#00A087FF", "#4DBBD5FF", "#E64B35FF", "#3C5488FF",
                    "#F38400", "#A1CAF1", "#BE0032", "#C2B280", "#848482", "#008856",
                    "#E68FAC", "#0067A5", "#604E97", "#F6A600", "#B3446C", "#DCD300",
                    "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26","#848482",
                    "#008856", "#E68FAC", "#0067A5", "#604E97", "#F6A600", "#B3446C",
                    "#DCD300", "#882D17")
  anno_color <- list()
  if (!is.null(anno_row)) {
    rowcolor <- color_palette[1:length(unique(anno_row[,1]))]
    names(rowcolor) <- unique(anno_row[,1])
    anno_color[[colnames(anno_row)]] <- rowcolor
  }
  if (!is.null(anno_col)) {
    colcolor <- color_palette[1:length(unique(anno_col[,1]))]
    names(colcolor) <- unique(anno_col[,1])
    anno_color[[colnames(anno_col)]] <- colcolor
  }
  pheatmap::pheatmap(ave, annotation_col = anno_col, annotation_row = anno_row, annotation_colors = anno_color)
}

#' multiheatmap function for showing heatmap of DEG
#' 
#' @examples 
#' diffg <- FindAllMarkers(pbmc, only.pos = T)
#' top20 <- diffg %>% group_by(cluster) %>% top_n(20, avg_logFC)
#' cell_order = unique(as.vector(top20$cluster))
#' pbmc <- ScaleData(pbmc, features = rownames(pbmc))
#' multiheatmap(object = pbmc, assay = 'RNA', DEGlist = top20, annotation = 'Annotation',
#'              group = 'Annotation', labelgene = NULL, labelgene_length = 2, slot = 'scale.data',
#'              mark_gene = F, block_gene = T, cell_order = unique(as.vector(top20$cluster)))
#' labelgene <- data.frame(gene = top20$gene, Group = top20$cluster, stringsAsFactors = F)
#' multiheatmap(object = pbmc, assay = 'RNA', DEGlist = top20, annotation = 'Annotation',
#'              group = 'Annotation', labelgene = labelgene, labelgene_length = 2, slot = 'scale.data',
#'              mark_gene = F, block_gene = T, cell_order = unique(as.vector(top20$cluster)))
#' multiheatmap(object = pbmc, assay = 'RNA', DEGlist = top20, annotation = 'Annotation',
#'              group = 'Annotation', labelgene = NULL, labelgene_length = 2, slot = 'scale.data',
#'              mark_gene = T, mark_genelist = c('CD4','CD8A','PPBP','CCR7','CD14','CD79A','NKG7'),
#'              cell_order = unique(as.vector(top20$cluster)))
#'
#' @export
multiheatmap <- function (object = NULL, assay = "RNA", DEGlist = NULL, annotation = NULL, 
    group = NULL, average.df = NULL, 
    cell_order = NULL, group_order = NULL, gene.order = NULL,
    slot = "scale.data", labelgene = NULL, labelgene_length = 2, 
    block_text_size = 10, mark_gene = F, mark_gene.fontsize = 10, 
    mark_genelist = NULL, cut.off = NULL,
    block_gene = T, block_alpha = 1,
    annotation_legend_name = "Maintypes", show_annotation_name = T, 
    annocolor = NULL, blockcolor = NULL, barcolor = NULL, legend.break = NULL,
    column_names_rot = 90, simple_anno_size = unit(0.5, "cm"), ...) 
{
    colnames(DEGlist)[which(colnames(DEGlist) == 'cluster')] <- "Annotation"
    if (all(sort(as.vector(unique(DEGlist$Annotation))) != 
            sort(as.vector(unique(object@meta.data[, annotation]))))) {
        stop("celltypes are different between DEGlist and object.")
    }
    a_g <- unique(object@meta.data[, c(annotation, group)])
    colnames(a_g) <- c("Annotation", "Group")
    a_g <- a_g[!duplicated(a_g$Annotation), ]
    DEGlist <- inner_join(DEGlist, a_g, by = "Annotation")
    DEGlist <- DEGlist[order(DEGlist$Group, DEGlist$Annotation), 
        ]
    if (!is.null(group_order)) {
        DEGlist$Group <- factor(DEGlist$Group, levels = group_order)
        DEGlist <- DEGlist[order(DEGlist$Group), ]
    }
    else {
        group_order <- unique(as.vector(DEGlist$Group))
    }
    if (is.null(cell_order)) {
        cell_order <- unique(as.vector(DEGlist$Annotation))
    }
    else {
        DEGlist$Annotation <- factor(DEGlist$Annotation, levels = cell_order)
        DEGlist <- DEGlist[order(DEGlist$Annotation), ]
    }
    DEGlist$Group <- factor(DEGlist$Group, levels = unique(DEGlist$Group))
    a_g$Annotation <- factor(a_g$Annotation, levels = cell_order)
    a_g <- a_g[order(a_g$Annotation), ]
    DEG <- unique(DEGlist[, c("Annotation", "gene")])
    if (!is.null(gene.order)) {
      paths <- gene.order
    }else{
      paths <- as.vector(DEG$gene)
    }
    missgene <- !paths %in% rownames(GetAssayData(object, slot = slot))
    if (any(missgene)) {
        stop(paste0("DEGs: ", paste(paths[missgene][1:10], collapse = ", "), 
            "... are not in the ", slot))
    }
    if (slot == "scale.data") {
        legend_title <- "Scaled expression"
    }
    else if (slot == "data") {
        legend_title <- "Normalized expression"
    }
    Idents(object) <- object@meta.data[, annotation]
    if (!is.null(average.df)) {
        data <- average.df
    }
    else {
        datalist <- AverageExpression(object, assays = assay, slot = slot)
        data <- datalist[[1]]
    }
    mat <- data[paths, , drop = FALSE]
    mat <- mat[, cell_order]
    rownames(mat) <- NULL
    mat <- as.matrix(mat)
    if (is.null(annocolor)) {
        num_g <- length(unique(object@meta.data[, group]))
        if (num_g < 8) {
            annocolor <- scPalette_anno_dark(num_g)
        }
        else {
            annocolor <- scPalette3(num_g)
        }
    }
    annodf <- data.frame(group = unique(as.vector(a_g$Group)), 
        color = annocolor, stringsAsFactors = F)
    colnames(annodf) <- c("Group", "color")
    annodf <- inner_join(a_g, annodf, by = "Group")
    rownames(annodf) <- annodf$Annotation
    annotationcol <- annodf$color
    names(annotationcol) <- annodf$Annotation
    ha = HeatmapAnnotation(annotation_group_name = cell_order, 
        show_legend = F, 
        name = annotation_legend_name, 
        show_annotation_name = show_annotation_name, 
        col = list(annotation_group_name = annotationcol),
        simple_anno_size = simple_anno_size)
    if (mark_gene) {
        pos <- which(paths %in% mark_genelist)
        row_anno = rowAnnotation(foo = anno_mark(at = pos, 
                                                 labels = as.vector(paths[pos]),
                                                 labels_gp = gpar(fontsize = mark_gene.fontsize))
                                 )
    }
    else if (block_gene) {
        if (is.null(labelgene)) {
            labelgene <- DEGlist %>% group_by(Annotation) %>% 
                top_n(3, avg_logFC)
            labelgene <- labelgene[, c("Annotation", "gene")]
            labelgene <- inner_join(labelgene, a_g, by = "Annotation")
        }
        labelgene_str <- sapply(unique(labelgene$Group), function(x) {
            genes <- as.vector(labelgene$gene[labelgene$Group == 
                x])
            genes <- unique(genes)
            num <- length(genes)
            s <- seq(from = 1, to = num, by = labelgene_length)
            s <- s[-1]
            genes[s] <- paste0("\n", genes[s])
            paste(genes, collapse = ", ")
        })
        if (is.null(blockcolor)) {
            block_num <- length(unique(labelgene$Group))
            if (block_num < 8) {
                blockcolor <- scPalette_anno_light(block_num)
            }
            else {
                blockcolor <- scPalette3(block_num)
            }
        }
        row_anno <- rowAnnotation(foo = anno_block(gp = gpar(fill = blockcolor, 
                                                             alpha = block_alpha), 
            labels = labelgene_str, labels_gp = gpar(col = "black", 
                fontsize = block_text_size), labels_rot = 0))
    }
    if (is.null(barcolor)) {
        barcolor <- c("#2596be","#93c4dc", "white", "#f5b08f", "#b93032", "#b11b2c")
        #barcolor <- c("blue", "white", "red")
    }
    if (!is.null(cut.off)) {
      if (cut.off[1] != -Inf) {
        mat[mat < cut.off[1]] <- cut.off[1]
        minvalue <- cut.off[1]
        minlabel <- paste0('<', cut.off[1])
      }else{
        minvalue <- floor(min(mat))
        minlabel <- floor(min(mat))
      }
      if (cut.off[2] != Inf) {
        mat[mat > cut.off[2]] <- cut.off[2]
        maxvalue <- cut.off[2]
        maxlabel <- paste0('>', cut.off[2])
      }else{
        maxvalue <- ceiling(max(mat))
        maxlabel <- ceiling(max(mat))
      }
    }
    if (!is.null(legend.break)) {
      if (!is.null(cut.off)) {
        legend.break <- c(minvalue, legend.break[2:(length(legend.break)-1)],maxvalue)
        labels <- c(minlabel,legend.break[2:(length(legend.break)-1)], maxlabel)
      }else{
        labels <- legend.break
        }
      col_fun = colorRamp2(legend.break, barcolor)
      hm <- Heatmap(mat, cluster_rows = F, cluster_columns = F, 
                    right_annotation = row_anno, top_annotation = ha, col = barcolor, 
                    row_split = DEGlist$Group, row_title = NULL, column_title = NULL, 
                    heatmap_legend_param = list(title = legend_title, 
                                                at = legend.break, 
                                                labels = labels,
                                                col_fun = col_fun), 
                    row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE, 
                    column_order = cell_order, column_names_rot = column_names_rot, ...)
    }else{
      hm <- Heatmap(mat, cluster_rows = F, cluster_columns = F, 
                    right_annotation = row_anno, top_annotation = ha, col = barcolor, 
                    row_split = DEGlist$Group, row_title = NULL, column_title = NULL, 
                    heatmap_legend_param = list(title = legend_title),
                    row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE, 
                    column_order = cell_order, column_names_rot = column_names_rot, ...)
    }
    return(hm)
}





#' Heatmap to show result of fish exact test of cell proportion.
#' 
#' @param testresult A data frame from result of cell_test.
#' @param label What type of text filled in heatmap when display_numbers is TRUE.
#' \itemize{
#'  \item{Pstar}: display '*' to show significance;
#'  \item{OR}: display values from 'OR' column.
#' }
#' 
#' @examples 
#' # 2x2
#' meta_epi <- meta200[meta200$maingroups == 'Epi_cells',]
#' statdf <- cell_test(meta_epi, Group = 'Group',Annotation = 'Annotation',
#' Sample = 'Sample', model = '2x2', test = 'fisher',alternative = 'two.sided',
#' adjust.method = 'BH')
#' # plot
#' plotOR(statdf, OR = 'OR', Status = 'Group', Celltype = 'Annotation', 
#' Pvalue = 'adj.Pvalue', display_numbers = F, label = 'Pstar',
#' log2.transfer = F, cut.off = c(-Inf,10), mid.point.positive = 1,coord.flip = F)

#' @export
plotOR <- function(testresult, OR = 'OR',
                   Status = 'Group', Celltype = 'Annotation', Pvalue = 'adj.Pvalue',
                   display_numbers = T, label = 'Pstar', log2.transfer = F,
                   negative.log10.transfer = F, 
                   rowname.order = NULL,
                   colname.order = NULL,
                   cut.off = NULL,  
                   breaks.hp = NA,
                   legend.breaks = NA, 
                   legend.breaks.raw = F,
                   legend.labels.raw = F,
                   colors = NULL, color.ramp = NULL,
                   mid.point.positive = 1, mid.point.negative = 0,
                   legend.title = NULL, coord.flip = F, ...) {
  all_result <- data.frame(Status = as.vector(testresult[,Status]),
                           Celltype = as.vector(testresult[,Celltype]),
                           OR = as.vector(testresult[,OR]),
                           Pvalue = as.vector(testresult[,Pvalue]),
                           stringsAsFactors = F)
  label_text <- reshape2::dcast(all_result[,c('Celltype','Status','Pvalue')],
                                formula = Celltype~Status, value.var = 'Pvalue')
  label_text <- data.frame(label_text, row.names = 1)
  label_star <- label_text
  label_star[label_text >= 0.05] <- ''
  label_star[label_text < 0.05] <- "\u2217"
  label_star[label_text < 0.01] <- "\u2217\u2217"
  label_star[label_text < 0.001] <- "\u2217\u2217\u2217"
  if (log2.transfer) {
    all_result$OR <- log2(all_result$OR + 0.00001)
  }else if (negative.log10.transfer) {
    all_result$OR <- -log10(all_result$OR + 0.00001)
  }
  OR <- reshape2::dcast(all_result[,c('Celltype','Status','OR')],
                        formula = Celltype~Status, value.var = 'OR')
  OR <- data.frame(OR, row.names = 1)
  if (!is.null(rowname.order)) {
    OR <- OR[rowname.order,]
    label_star <- label_star[rowname.order,]
  }
  if (!is.null(colname.order)) {
    OR <- OR[,colname.order]
    label_star <- label_star[,colname.order]
  }
  if (!is.null(cut.off)) {
    OR[OR < cut.off[1]] <-  cut.off[1]
    OR[OR > cut.off[2]] <-  cut.off[2]
  }
  paletteLength <- 100
  if (is.null(colors)) {
    colors <- (grDevices::colorRampPalette(rev(c( "#D73027", "#F46D43",'white', "#74ADD1", '#4575B4'))))(paletteLength)
  }else{
    colors <- grDevices::colorRampPalette(colors)(paletteLength)
  }
  if (!is.null(color.ramp)) colors <- color.ramp
  rangel <- c(floor(min(OR)), floor(max(OR)))
  if (any(c(-Inf, Inf) %in% rangel)) {
    stop('-Inf/Inf was found in plotdata. Please set cut off.')
  }
  if (mid.point.positive < rangel[1] | mid.point.positive > rangel[2]) {
    message('Warning: mid.point.positive was out of range, replaced by average value of OR.')
    mid.point.positive <- mean(rangel)
  }
  if (rangel[1] < 0) {
    length.all <- 7
    step.length <- ceiling((rangel[2] - rangel[1])/length.all)
    legend_breaks <- c(rev(-seq(mid.point.negative, abs(rangel[1]), step.length)[-1]),
                       seq(mid.point.negative, rangel[2], step.length), max(OR))
    breaks <- c(seq(min(OR), mid.point.negative, length.out=ceiling(paletteLength/2) + 1),
                seq(max(OR)/paletteLength, max(OR), length.out=floor(paletteLength/2)))
    if (is.null(legend.title)) {
      legend.title <- 'log2(OR)'
    }
  }else{
    length.all <- 6
    
    if (rangel[2] <= 1) {
      legend_breaks <- c(seq(rangel[1],max(OR), by = 0.5),max(OR))
      breaks <- c(seq(0, mid.point.positive, length.out = ceiling(paletteLength/2)),
                  seq(mid.point.positive, max(OR), length.out=floor(paletteLength/2))[-1])
    }else{
      step.length <- ceiling((rangel[2] -rangel[1])/length.all)
      legend_breaks <- c(0,mid.point.positive, round(seq(mid.point.positive, rangel[2], step.length),0), max(OR))
      breaks <- c(seq(0, mid.point.positive, length.out = ceiling(paletteLength/2)),
                  seq(mid.point.positive, max(OR), length.out=floor(paletteLength/2))[-1])
    }
    if (is.null(legend.title)) {
      legend.title <- 'Odds ratio'
    }
  }
  if (any(!is.na(legend.breaks))) {legend_breaks <- legend.breaks}
  len <- length(legend_breaks)
  if (!is.null(cut.off)) {
    if (cut.off[1] == -Inf) {
      min.legend.break <- legend_breaks[1]
    }else{
      min.legend.break <- paste0('<',legend_breaks[1])
    }
    if (cut.off[2] == Inf) {
      max.legend.break <- legend_breaks[len-1]
    }else{
      max.legend.break <- paste0('<',legend_breaks[len-1])
    }
    legend_labels <- c(min.legend.break,legend_breaks[-c(1,len-1, len)],
                       max.legend.break, paste0(legend.title, "\n"))
  }else{
    legend_labels <- c(legend_breaks[-len],paste0(legend.title, "\n"))
  }
  if (legend.labels.raw) {legend_labels <- NA}
  if (legend.breaks.raw) {legend_breaks <- NA}
  if (any(!is.na(breaks.hp))) {breaks <- breaks.hp}
  #print(legend_breaks)
  #print(legend_labels)
  if (coord.flip) {
    OR <- as.data.frame(t(OR))
    label_star <- as.data.frame(t(label_star))
  }
  if (display_numbers) {
    if (label == 'Pstar') {
      display_numbers = label_star
    }else if (label == 'OR') {
      display_numbers = T
    }
  }
  breaks <- unique(breaks)
  pheatmap::pheatmap(OR, display_numbers = display_numbers,
                     treeheight_row = 10, treeheight_col = 10,
                     color = colors,
                     breaks = breaks,
                     legend_breaks = legend_breaks,
                     legend_labels = legend_labels,...)
}
#' Heatmap of result from MASC.
#' 
#' @examples 
#' stat2x2 <- MASC_multi(meta200, cluster = 'maingroups', contrast = 'Group',
#' random_effects = 'Sample',model = '2x2')
#' plotMASC(stat2x2)
#' 
#' @export
plotMASC <- function (all_result, OR = "OR", display_numbers = T, label = "Pvalue", 
                      log2.transfer = F, legend.title = NULL, coord.flip = F) 
{
  label_text <- reshape2::dcast(all_result[, c("Celltype", "Status", "Pvalue")],
                                formula = Celltype ~ Status)
  label_text <- data.frame(label_text, row.names = 1)
  label_star <- label_text
  label_star[label_text >= 0.05] <- ""
  label_star[label_text < 0.05] <- "∗"
  label_star[label_text < 0.01] <- "∗∗"
  label_star[label_text < 0.001] <- "∗∗∗"
  if (log2.transfer) {
    all_result$OR <- log2(all_result$OR)
  }
  OR <- reshape2::dcast(all_result[, c("Celltype", "Status", 
                                       "OR")], formula = Celltype ~ Status)
  OR <- data.frame(OR, row.names = 1)
  paletteLength <- 100
  colors <- (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9,name = "RdYlBu"))))(paletteLength)
  colors <- (grDevices::colorRampPalette(rev(c("#D73027", "#F46D43", 
                                               "white", "#74ADD1", "#4575B4"))))(100)
  rangel <- c(floor(min(OR)), floor(max(OR)))
  if (rangel[1] < 0) {
    length.all <- 7
    step.length <- ceiling((rangel[2] - rangel[1])/length.all)
    legend_breaks <- c(rev(-seq(0, abs(rangel[1]), step.length)[-1]), 
                       seq(0, rangel[2], step.length), max(OR))
    breaks <- c(seq(min(OR), 0, length.out = ceiling(paletteLength/2) + 
                      1), seq(max(OR)/paletteLength, max(OR), length.out = floor(paletteLength/2)))
    if (is.null(legend.title)) {
      legend.title <- "log2(OR)"
    }
  }
  else {
    length.all <- 6
    step.length <- ceiling((rangel[2] - rangel[1])/length.all)
    legend_breaks <- c(0, 1, round(seq(1, rangel[2], step.length), 
                                   0), max(OR))
    breaks <- c(seq(0, 1, length.out = ceiling(paletteLength/2)), 
                seq(1, max(OR), length.out = floor(paletteLength/2))[-1])
    if (is.null(legend.title)) {
      legend.title <- "Odds ratio"
    }
  }
  legend_labels <- c(legend_breaks[-length(legend_breaks)], 
                     paste0(legend.title, "\n"))
  if (coord.flip) {
    OR <- as.data.frame(t(OR))
    label_star <- as.data.frame(t(label_star))
  }
  if (display_numbers) {
    if (label == "Pvalue") {
      display_numbers = label_star
    }
    else if (label == "OR") {
      display_numbers = T
    }
  }
  pheatmap::pheatmap(OR, display_numbers = display_numbers, 
                     treeheight_row = 10, treeheight_col = 10, color = colors, 
                     breaks = breaks, legend_breaks = legend_breaks,
                     legend_labels = legend_labels, ...)
}
#' Scale color of heatmap for pheatmap.
#' 
#' @param OR matrix of value. The format is features x samples
#' @param mid.point.positive Value to indicate the middle point of legend. 
#' And all the value in matrix are positive (>0).
#' @param mid.point.negative Value to indicate the middle point of legend. 
#' And the minimum value in matrix is negative (<0).
#' 
#' @examples 
#' 
#' @export
scale_color_heatmap <- function(OR, 
                                colors = NULL, 
                                mid.point.positive = 1, mid.point.negative = 0, 
                                length.node.num = NULL,
                                legend.breaks = NULL, cut.off = NULL,
                                legend.title = NULL, paletteLength = 100) {
  
  if (is.null(colors)) {
    colors <- grDevices::colorRampPalette(rev(c("#D73027", 
                                                "#F46D43", "white", "#74ADD1", "#4575B4")))(paletteLength)
  }
  else {
    colors <- (grDevices::colorRampPalette(colors))(paletteLength)
  }
  rangel <- c(floor(min(OR)), floor(max(OR)))
  if (any(c(-Inf, Inf) %in% rangel)) {
    stop("-Inf/Inf was found in plotdata. Please set cut off.")
  }
  if (mid.point.positive < rangel[1] | mid.point.positive > 
      rangel[2]) {
    message("Warning: mid.point.positive was out of range, replaced by average value of OR.")
    mid.point.positive <- mean(rangel)
  }
  if (rangel[1] < 0) {
    if (is.null(length.node.num)) length.node.num <- 7
    step.length <- ceiling((rangel[2] - rangel[1])/length.node.num)
    legend_breaks <- c(rev(-seq(-mid.point.negative, abs(rangel[1]), step.length)[-1]),
                       seq(mid.point.negative, rangel[2], step.length), max(OR))
    #breaks <- c(seq(min(OR), mid.point.negative, length.out = ceiling(paletteLength/2) + 1),
    #            seq(max(OR)/paletteLength, max(OR), length.out = floor(paletteLength/2)))
    breaks <- c(seq(min(OR), mid.point.negative, length.out = ceiling(paletteLength/2) + 1),
                seq(mid.point.negative, max(OR), length.out = floor(paletteLength/2))[-1])
    if (is.null(legend.title)) {
      legend.title <- "log2(OR)"
    }
  }
  else {
    if (is.null(length.node.num)) length.node.num <- 6
    if (rangel[2] <= 1) {
      legend_breaks <- c(seq(rangel[1], max(OR), by = 0.5), 
                         max(OR))
      breaks <- c(seq(0, mid.point.positive, length.out = ceiling(paletteLength/2)), 
                  seq(mid.point.positive, max(OR), length.out = floor(paletteLength/2))[-1])
    }
    else {
      step.length <- ceiling((rangel[2] - rangel[1])/length.node.num)
      legend_breaks <- c(0, mid.point.positive, round(seq(mid.point.positive, 
                                                          rangel[2], step.length), 0), max(OR))
      breaks <- c(seq(0, mid.point.positive, length.out = ceiling(paletteLength/2)), 
                  seq(mid.point.positive, max(OR), length.out = floor(paletteLength/2))[-1])
    }
    if (is.null(legend.title)) {
      legend.title <- "Odds ratio"
    }
  }
  if (!is.null(legend.breaks)) {
    legend_breaks <- legend.breaks
  }
  len <- length(legend_breaks)
  if (!is.null(cut.off)) {
    if (cut.off[1] == -Inf) {
      min.legend.break <- legend_breaks[1]
    }
    else {
      min.legend.break <- paste0("<", legend_breaks[1])
    }
    if (cut.off[2] == Inf) {
      max.legend.break <- legend_breaks[len - 1]
    }
    else {
      max.legend.break <- paste0("<", legend_breaks[len - 
                                                      1])
    }
    legend_labels <- c(min.legend.break, legend_breaks[-c(1, len - 1, len)], 
                       max.legend.break, paste0(legend.title,"\n"))
  }
  else {
    legend_labels <- c(legend_breaks[-len], paste0(legend.title,  "\n"))
  }
  legendlist <- list(color = colors, 
                     breaks = breaks,
                     legend_breaks = legend_breaks, 
                     legend_labels = legend_labels)
  return(legendlist)
}
