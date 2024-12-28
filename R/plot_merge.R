#' Scater plot of two genes.
#'
#' @param SeuratS4 a Seurat object.
#' @param cl select which cluster to plot.
#' @param gene1 one of two genes.
#' @param gene2 one of two genes.
#' @examples
#' data('pbmc')
#' plotcogene(pbmc,'Naive CD4 T','CD4','CCR7')
#' @export
plotcogene <- function(SeuratS4, cl, gene1, gene2, return.frac = F) {
  data <- as.data.frame(t(as.matrix(SeuratS4@assays$RNA@data[c(gene1, gene2),
                                                             Idents(SeuratS4) == cl])))
  p <- ggplot(data, aes(data[, gene1], data[, gene2])) + geom_point() +
    xlab(gene1) + ylab(gene2) + labs(title = paste0("Cluster ", cl))
  frac <- length(which(apply(data, 1, function(x) {
    length(which(x > 0)) == 2
  })))/nrow(data)
  print(paste0("Co-expression fraction: ", frac))
  if (return.frac) {
    return(frac)
  }else{
    return(p)
  }
}

#'Gradual color plot of two co-expression genes.
#'
#' @param exprmat a normalize 2 selected genes expression dataframe(cell x gene).
#' @param xy cell location from seurat_object@@reduction$tsne@@cell.embeddings.
#' @examples
#' data('pbmc')
#' marker_to_plot <- c("CD4","CCR7")
#' exp <- log2(as.data.frame(t(pbmc@@assays$RNA@@data[marker_to_plot,]))+1)
#' xy <- pbmc@@reductions$umap@@cell.embeddings
#' mergeplot(xy,exp)
#' @export

mergeplot <- function(xy, exprmat) {
  if (ncol(exprmat) == 2) {
    # exprmat <- log(exprmat + 1)
    RGmat <- apply(exprmat, 2, function(x) (x - min(x))/(max(x) - min(x)))
    RGBmat <- cbind(RGmat, 0)
    color <- apply(RGBmat, 1, function(x) rgb(x[1], x[2], x[3]))
    # color <- gsub('#FFFF00','#DCDCDC',color) arrangesum is used for
    # arrange data
    arrangesum <- apply(RGBmat, 1, sum)
    xymat <- data.frame(xy, color, arrangesum)
    xynames <- colnames(xy)
    colnames(xymat) <- c(xynames, "RGBcolor", "arr")
    # plot.data <- plyr::arrange(xymat,arr)
    plot.data <- xymat[order(xymat$arr, decreasing = F), ]
    p <- ggplot(plot.data, aes(x = plot.data[, 1], y = plot.data[,
                                                                 2], colour = RGBcolor)) + geom_point(size = 0.1, color = plot.data$RGBcolor) +
      theme_bw() + theme(panel.grid.major = element_line(colour = NA),
                         panel.grid.minor = element_line(colour = NA)) + theme(legend.position = "none")
    print("Plot genes: ")
    print(colnames(exprmat))
    print(p)
    ## generate color patttern
    polar2 <- apply(exprmat, 2, function(x) {
      seq(min(x), max(x), length.out = 51)
    })
    color51 <- seq(0, 255, 5)
    rowlist <- list()
    for (i in 1:51) {
      for (j in 1:51) {
        color <- rgb(color51[i], color51[j], 0, maxColorValue = 255)
        rowlist[[paste0(i, "_", j)]] <- c(polar2[i, 1], polar2[j,
                                                               2], color)
      }
    }
    colorfill <- data.frame(rowlist)
    colorfill <- as.data.frame(t(colorfill))
    colorfill[, c(1, 2)] <- lapply(colorfill[, c(1, 2)], function(x) {
      as.numeric(as.character(x))
    })
    sapply(colorfill, class)
    colnames(colorfill) <- c(colnames(exprmat), "RGBcolor")
    p <- ggplot(colorfill, aes(x = colorfill[, 1], y = colorfill[,
                                                                 2])) + geom_point(size = 0.5, color = colorfill$RGBcolor) +
      theme_bw() + theme(panel.grid.major = element_line(colour = NA),
                         panel.grid.minor = element_line(colour = NA)) + theme(legend.position = "none")
    print(p)
  } else {
    print("Only TWO genes allowed.")
  }
}

#'3 color plot of two co-expression genes.
#'
#' @param SeuratS4 a Seurat object.
#' @param marker_to_plot two Genes to plot.
#' @examples
#' data('pbmc')
#' marker_to_plot <- c("CD4","CCR7")
#' mergeplot1(pbmc,marker_to_plot)
#' @export

## 3 fill bar
mergeplot1 <- function(SeuratS4, marker_to_plot, active.assay = 'RNA', reduction = 'umap') {
  xy <- SeuratS4@reductions[[reduction]]@cell.embeddings
  if (length(marker_to_plot) == 2) {
    if (sum(marker_to_plot %in% rownames(SeuratS4@assays[[active.assay]]@data)) !=
        2) {
      print("Input gene doesn't exist in data!\nExist.")
      break
    } else {
      exprmat <- log2(as.data.frame(t(as.matrix(SeuratS4@assays$RNA@data[marker_to_plot,
      ]))) + 1)
      RGmat <- apply(exprmat, 2, function(x) (x - min(x))/(max(x) -
                                                             min(x)))
      # color <- rep('#D3D3D3',nrow(RGmat))
      color <- apply(RGmat, 1, function(x) {
        if (x[1] > 0 & x[2] == 0) {
          rgb(255, 200 * (1 - x[1]), 200 * (1 - x[1]), maxColorValue = 255)
        } else if (x[1] == 0 & x[2] > 0) {
          rgb(200 * (1 - x[2]), 255, 200 * (1 - x[2]), maxColorValue = 255)
        } else if (x[1] == 0 & x[2] == 0) {
          # '#D3D3D3'
          rgb(211, 211, 211, maxColorValue = 255)
        } else {
          # yellow
          rgb(255, 255, 200 * (1 - min(x)), maxColorValue = 255)
          # blue rgb(200*(1-min(x)),200*(1-min(x)),255,maxColorValue = 255)
        }
      })
      comm_pos <- apply(RGmat, 1, function(x) x[1] > 0 & x[2] > 0)
      comm_min <- apply(RGmat[comm_pos, ], 1, min)
      # comm_min <- log(comm_min+1)
      scalec <- (comm_min - min(comm_min))/(max(comm_min) - min(comm_min))
      comm_col <- c()
      for (x in scalec) {
        ## blue comm_col <- c(comm_col,rgb(200*(1-x),200*(1-x),255,maxColorValue
        ## = 255)) yellow
        comm_col <- c(comm_col, rgb(255, 255, 200 * (1 - min(x)),
                                    maxColorValue = 255))
      }
      color[comm_pos] <- comm_col
      ## arrangesum is used for arrange data
      arrangesum <- apply(RGmat, 1, sum)
      xymat <- data.frame(xy, color, arrangesum)
      xynames <- colnames(xy)
      colnames(xymat) <- c(xynames, "RGBcolor", "arr")
      plot.data <- plyr::arrange(xymat, arr)
      # plot.data <- xymat[order(xymat$arr,decreasing = F),]
      p <- ggplot(plot.data, aes(x = plot.data[, 1], y = plot.data[,
                                                                   2], color = RGBcolor)) + geom_point(size = 0.5, color = plot.data$RGBcolor) +
        theme_bw() + theme(panel.grid.major = element_line(colour = NA),
                           panel.grid.minor = element_line(colour = NA)) + # theme(legend.position = 'none')+
        xlab(xynames[1]) + ylab(xynames[2])
      print("Plot genes: ")
      print(colnames(exprmat))
      print(p)
    }
  } else {
    print("Only TWO genes allowed.")
  }
}
#'3 color plot of two co-expression genes.
#'
#' @param SeuratS4 a Seurat object.
#' @param marker_to_plot two Genes to plot.
#' @param expr.cut.off Expression below Cut off will be assigned to 0.
#' @examples
#' data('pbmc')
#' marker_to_plot <- c("CD4","LEF1")
#' mergeplot2(pbmc,marker_to_plot, expr.cut.off = 0.1)
#' @export
mergeplot2 <- function(SeuratS4, marker_to_plot, active.assay = 'RNA',
                       reduction = 'umap', expr.cut.off = NULL) {
  xy <- SeuratS4@reductions[[reduction]]@cell.embeddings
  if (length(marker_to_plot) == 2) {
    if (sum(marker_to_plot %in% rownames(SeuratS4@assays[[active.assay]]@data)) !=
        2) {
      print("Input gene doesn't exist in data!\nExist.")
      break
    } else {
      mat <- SeuratS4@assays$RNA@data[marker_to_plot,]
      if (!is.null(expr.cut.off)) {
        mat[mat < expr.cut.off] <- 0
      }
      exprmat <- log2(as.data.frame(t(as.matrix(mat))) + 1)
      RGmat <- apply(exprmat, 2, function(x) (x - min(x))/(max(x) -
                                                             min(x)))
      # color <- rep('#D3D3D3',nrow(RGmat))
      color <- apply(RGmat, 1, function(x) {
        if (x[1] > 0 & x[2] == 0) {
          rgb(255, 200 * (1 - x[1]), 200 * (1 - x[1]), maxColorValue = 255)
        } else if (x[1] == 0 & x[2] > 0) {
          rgb(200 * (1 - x[2]), 255, 200 * (1 - x[2]), maxColorValue = 255)
        } else if (x[1] == 0 & x[2] == 0) {
          # '#D3D3D3'
          rgb(211, 211, 211, maxColorValue = 255)
        } else {
          # yellow
          rgb(255, 255, 200 * (1 - min(x)), maxColorValue = 255)
          # blue rgb(200*(1-min(x)),200*(1-min(x)),255,maxColorValue = 255)
        }
      })
      comm_pos <- apply(RGmat, 1, function(x) x[1] > 0 & x[2] > 0)
      comm_min <- apply(RGmat[comm_pos, ], 1, min)
      # comm_min <- log(comm_min+1)
      scalec <- (comm_min - min(comm_min))/(max(comm_min) - min(comm_min))
      comm_col <- c()
      for (x in scalec) {
        ## blue comm_col <- c(comm_col,rgb(200*(1-x),200*(1-x),255,maxColorValue
        ## = 255)) yellow
        comm_col <- c(comm_col, rgb(255, 255, 200 * (1 - min(x)),
                                    maxColorValue = 255))
      }
      color[comm_pos] <- comm_col
      ## arrangesum is used for arrange data
      arrangesum <- apply(RGmat, 1, sum)
      xymat <- data.frame(xy, color, arrangesum)
      xynames <- colnames(xy)
      colnames(xymat) <- c(xynames, "RGBcolor", "arr")
      plot.data <- plyr::arrange(xymat, arr)
      # plot.data <- xymat[order(xymat$arr,decreasing = F),]
      p <- ggplot(plot.data, aes(x = plot.data[, 1], y = plot.data[,
                                                                   2], color = RGBcolor)) + geom_point(size = 0.5, color = plot.data$RGBcolor) +
        theme_bw() + theme(panel.grid.major = element_line(colour = NA),
                           panel.grid.minor = element_line(colour = NA)) + # theme(legend.position = 'none')+
        xlab(xynames[1]) + ylab(xynames[2])
      print("Plot genes: ")
      print(colnames(exprmat))
      print(p)
    }
  } else {
    print("Only TWO genes allowed.")
  }
}
