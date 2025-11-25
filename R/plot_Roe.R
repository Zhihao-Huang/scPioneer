#' Plot Ro/e heatmap
#' @author Li Ziyi
#' @export
ROIE_plot <- function(ROIE_table, max = 2.5, font_size = 12, geom_text_size = 11){
  ROIE_table <- melt(ROIE_table)
  ROIE_table$text <- round(ROIE_table$value,2)
  ROIE_table$value[ROIE_table$value > max] <- max
  #ROIE_table$value[ROIE_table$value == 0] <- NA
  #ROIE_table$text[ROIE_table$text == 0] <- "-"
  p <- 
    ggplot(data = ROIE_table, aes(Var2, Var1, fill = value)) +
    geom_tile() + 
    scale_fill_gradientn(name = "Ro/e", colours = colorRampPalette(brewer.pal(4, "YlOrRd"))(100), na.value = "lightgrey") +
    geom_text(aes(label = text), size = geom_text_size  / .pt) +
    labs(x = "", y = "") +
    theme_cowplot(font_size = font_size) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
          axis.line = element_blank(), 
          axis.ticks = element_blank()) +
    guides(fill = guide_colourbar(barwidth = 1, barheight = 10))
  return(p)
}
#' Plot Ro/e heatmap
#' @export
plot_Roe <- function(meta, colname, rowname, cluster_col = F, cluster_row = T,
                     max.Roe = 2.5, font_size = 12,
                     x.size = 10, y.size = 10, geom_text_size = 11) {
  rovemat <- rove(x=meta[, rowname],
                  y=meta[, colname],plot=F)
  data <- pheatmap::pheatmap(rovemat, silent = T)
  
  if (cluster_col) {
    g <- data$tree_col$labels[data$tree_col$order]
    rovemat <- rovemat[,g]
  }
  if (cluster_row) {
    cell.r <- data$tree_row$labels[data$tree_row$order]
    rovemat <- rovemat[cell.r,]
  }
  p3_1 <- ROIE_plot(rovemat, max = max.Roe, font_size = font_size,
                    geom_text_size = geom_text_size) + 
    theme(axis.text.x  = element_text(size = x.size),
          axis.text.y  = element_text(size = y.size))
  p3_1$data$Var1 <- factor(p3_1$data$Var1, levels = rev(levels(p3_1$data$Var1)))
  return(p3_1)
}

#' Heatmap to show result of Ro/e matrix.
#' 
#' @param label What type of text filled in heatmap when display_numbers is TRUE.
#' \itemize{
#'  \item{Pstar}: display '-/+' to show significance;
#'  \item{OR}: display values from 'OR' column.
#' }
#' 
#' @examples 
#' # plot
#' plotOR(statdf, OR = 'OR', Status = 'Group', Celltype = 'Annotation', 
#' Pvalue = 'adj.Pvalue', display_numbers = F, label = 'Pstar',
#' log2.transfer = F, cut.off = c(-Inf,10), mid.point.positive = 1,coord.flip = F)

#' @export
plotRO <- function(rovemat,
                   display_numbers = T, label = 'Pstar', log2.transfer = F,
                   negative.log10.transfer = F, 
                   paletteLength = 100,
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
  label_text <- rovemat
  label_star <- label_text
  #label_star[label_text >= 0] <- '±'
  label_star[label_text >= 0] <- "±"
  label_star[label_text > 1] <- "+"
  label_star[label_text > 1.5] <- "++"
  label_star[label_text > 3] <- "+++"
  if (log2.transfer) {
    label_text <- log2(label_text + 0.00001)
  }else if (negative.log10.transfer) {
    label_text <- -log10(label_text + 0.00001)
  }
  OR <- label_text
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
      legend.title <- 'Ro/e'
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
