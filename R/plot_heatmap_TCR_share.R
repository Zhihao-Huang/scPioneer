shareTCRdata_order_old <- function(mat) {
  typenum <- ncol(mat)
  #celltype <- gsub('^.*_','',colnames(mat))
  celltype <- colnames(mat)
  if (typenum == 2) {
    data <- mat[order(mat[,1],mat[,2],decreasing = T),]
    legenddata <- data.frame(shareTCR = paste0(celltype[1], ' and ', celltype[2]),
                             y = c(nrow(mat)),
                             stringsAsFactors = F)
  }
  if (typenum == 3) {
    #order in 3 segments.
    mat1 <- mat[mat[,3] == 0, ,drop = F]
    mat2 <- mat[mat[,1] == 0, ,drop = F]
    mat3 <- mat[mat[,2] == 0, ,drop = F]
    mat4 <- mat[apply(mat, 1, function(x) sum(x > 0) == 3),]
    #order
    mat1 <- mat1[order(mat1[,1],decreasing = T),,drop = F]
    mat2 <- mat2[order(mat2[,2],decreasing = T),,drop = F]
    mat3 <- mat3[order(mat3[,3],decreasing = T),,drop = F]
    if (dim(mat4)[1] != 0) {
      mat4 <- mat4[order(mat4[,1],mat4[,2],mat4[,3],decreasing = T),]
    }
    data <- rbind(mat1,mat2,mat3,mat4)
    legenddata <- data.frame(shareTCR = c(paste0(celltype[1], ' and ', celltype[2]),
                                          paste0(celltype[2], ' and ', celltype[3]),
                                          paste0(celltype[1], ' and ', celltype[3]),
                                          paste0(celltype[1],', ',celltype[2] , ' and ', celltype[3])),
                             y = c(nrow(mat1),nrow(mat2),nrow(mat3),nrow(mat4)))
  }
  if (typenum == 4) {
    #order pairs with 2 types.
    #mat
    matlist <- list()
    #legenddata
    shareTCR <- c()
    y <- c()
    #order
    for (i in 1:(typenum - 1)) {
      for (j in (i + 1) : typenum) {
        othertype <- !1:typenum %in% c(i,j)
        zeropos <- apply(mat[,othertype], 1, function(x) sum(x == 0) == sum(othertype))
        mat_ij <- mat[zeropos, ,drop = F]
        mat_ij <- mat_ij[order(mat_ij[,i],decreasing = T),,drop = F]
        pairsname <- paste0(colnames(mat)[i] , ' and ', colnames(mat)[j])
        if (dim(mat_ij)[1] != 0) {
          matlist[[pairsname]] <- mat_ij
          shareTCR <- c(shareTCR, pairsname)
          y <- c(y, nrow(mat_ij))
        }
      }
    }
    #order pairs with 3 types.
    for (i in rev(1:typenum)) {
      mat1 <- mat[mat[,i] == 0 & apply(mat, 1, function(x) sum(x == 0) == 1), ,drop=F]
      if (dim(mat1)[1] != 0) {
        pos <- which(! 1:typenum %in% i)
        mat1 <- mat1[order(mat1[,pos[1]],decreasing = T),,drop = F]
        if (dim(mat1)[1] != 0) {
          pairsname <- paste0(celltype[pos[1]] , ', ',celltype[pos[2]], ' and ', celltype[pos[3]])
          matlist[[pairsname]] <- mat1
          shareTCR <- c(shareTCR, pairsname)
          y <- c(y, nrow(mat1))
        }
      }
    }
    #order pairs with 4 types.
    mat4 <- mat[apply(mat, 1, function(x) sum(x > 0) == 4),,drop = F]
    if (dim(mat4)[1]!=0) {
      mat4 <- mat4[order(mat4[,1],mat4[,2],mat4[,3],mat4[,4],decreasing = T),,drop = F]
      pairsname <- 'All types'
      matlist[[pairsname]] <- mat4
      shareTCR <- c(shareTCR, pairsname)
      y <- c(y, nrow(mat4))
    }
    data <- do.call(rbind, matlist)
    legenddata <- data.frame(shareTCR = shareTCR,
                             y = y,
                             stringsAsFactors = F)
  }
  return(list(data = data, legdata = legenddata))
}
#' order heatmap data and get share number of TCR for length of annotation bar.
shareTCRdata_order <- function(mat) {
  df <- as.data.frame.matrix(mat)
  # data
  matlist <- list()
  # legenddata
  celltype <- colnames(df)
  shareTCR <- c()
  y <- c()
  # choose all tuples of column index
  cell_num <- ncol(df)
  tuples <- unlist(lapply(1:cell_num, function(n) combn(1:cell_num, n, simplify=FALSE)), recursive = F)
  # remove individual one
  tuples <- tuples[-c(1:cell_num)]
  # get all possible groups of mats
  for (x in tuples) {
    groups <- colnames(df)[x]
    groups_0 <- colnames(df)[!colnames(df) %in% groups]
    pos <- apply(df[, groups, drop = F],1, function(x) !any(x == 0)) &
      apply(df[, groups_0, drop = F],1, function(x) !any(x > 0))
    if (!any(pos)) {next}
    subdf <- df[pos,]
    subdf <- subdf %>% arrange(desc(subdf))
    anno_text <- paste(groups, collapse = ',')
    shareTCR <- c(shareTCR, anno_text)
    y <- c(y, nrow(subdf))
    #### 20220822 add clone name as row names. 
    #### 20220822.01 stat
    rownames(subdf) <- paste0(anno_text, '.', rownames(subdf))
    #### 20220822.01 end
    matlist[[anno_text]] <- subdf
  }
  #### 20220822.02 stat
  names(matlist) <- NULL
  #### 20220822.02 end
  orderdf <- do.call(rbind, matlist)
  legenddata <- data.frame(shareTCR = shareTCR,
                           y = y,
                           stringsAsFactors = F)
  return(list(data = orderdf, legdata = legenddata))
}
#' Heatmap of share TCR
#' 
#' @examples 
#' T_hgsoc <- meta200
#' typemat1 <- data.frame(Tissue = c('PT','PT'),
#' celltype = c('T08','T10'))
#' list1 <- shareTCR_heatmap(meta = T_hgsoc, group = 'Group_abb', annotation = 'Annotation_abb',
#' clonotype_index = 'clonotype_index', typemat = typemat1, return.data = T)
#' list1$plot
#' 
#' @export
shareTCR_heatmap <- function(meta, group, annotation, clonotype_index, typemat,
                             cell.order = NULL,
                             keep.group = NULL,
                             show.unique.TCR.group = NULL,
                             legend.ncol = NULL,
                             widths = c(1,8,6), x.size = 15,
                             legend.heatmap.size = 10, legend.anno.size = 10,
                             legend.heatmap.title.size = 15,
                             legend.anno.title.size = 15,
                             return.data = F) {
  colnames(typemat) <- c('Tissue','celltype')
  #if (!nrow(typemat) %in% c(2,3,4)) {
  #  stop(paste0('Only 2~4 groups were allowed, but ',nrow(typemat),' groups were provided.'))
  #}
  if(any(!typemat$Tissue %in% meta[,group]) & any(unique(typemat$Tissue) != 'all')){
    stop("Some tissues in typemat didn't exist in meta group.")
  }
  if(any(!typemat$celltype %in% meta[,annotation])) {
    stop("Some celltypes in typemat didn't exist in meta annotation.")
  }
  # use all tissues.
  if (!any(unique(typemat$Tissue) != 'all')) {
    meta$index <- meta[,annotation]
    select_g <- typemat$celltype
  }else{
    meta$index <- paste0(meta[,group],'_',meta[,annotation])
    select_g <- paste0(typemat$Tissue, '_',typemat$celltype)
  }
  #select
  meta <- meta[meta$index %in% select_g,]
  metac <- meta[,c(group, clonotype_index,'index')]
  metac <- metac[!is.na(metac[,clonotype_index]),]
  mat <- as.matrix(table(metac[,c(clonotype_index, 'index')]))
  
  # only keep TCR from one group.(remove empty part)
  if (!is.null(keep.group)) {
    df <- as.data.frame.matrix(mat)
    pos <- df[,keep.group] == 0
    mat <- mat[!pos, ]
  }
  #only share TCR
  sharepos <- apply(mat,1,function(x) sum(x>0) > 1)
  mat2 <- mat[sharepos,]
  if (sum(sharepos) == 0) {
    message('No shared TCR found.')
    return(list(stat = NA, plot = NA))  
  }
  #order mat and get legend data;
  mat_leglist <- shareTCRdata_order(mat2)
  legdata <- mat_leglist[['legdata']]
  data <- mat_leglist[['data']]
   
  # show all TCR of one type (including unique TCR)
  if (!is.null(show.unique.TCR.group)) {
    # add unique counts to data
    df <- as.data.frame.matrix(mat)
    uniq.pos <- apply(df[,!colnames(df) %in% c(show.unique.TCR.group)], 1, function(x) sum(x) == 0)
    df <- df[uniq.pos,]
    rownames(df) <- paste0(paste(show.unique.TCR.group,collapse = ','),'.',
                           rownames(df))
    df <- df[order(df[,show.unique.TCR.group, drop = F], decreasing = T),]
    data <- rbind(data, df)
    # add unique anno to legdata
    annodf <- data.frame(shareTCR = paste(show.unique.TCR.group,collapse = ','),
                         y = nrow(df), # y-axis of annotation bar
                         stringsAsFactors = F)
    # combine
    legdata <- rbind(legdata, annodf)
  }
  datalist <- list()
  datalist[['stat']] <- data
  # transform to data frame, matm, for ploting purpose.
  data[data>8] <- 9
  matm <- data.frame(colonalID = rep(rownames(data),ncol(data)),
                     T_celltype = rep(colnames(data),each = nrow(data)),
                     Cells = unlist(as.data.frame(data)))
  matm$celltype <- matm$T_celltype
  if (is.null(cell.order)) {
    cell.order <- unique(matm$celltype)
  }
  matm$celltype <- factor(matm$celltype, levels = cell.order)
  matm$Cells <- as.character(matm$Cells)
  if (any(matm$Cells > 8)) {
    matm$Cells <- gsub('9','>8',matm$Cells)
    matm$Cells <- factor(matm$Cells, levels = rev(c(as.character(0:8),'>8')))
  }else{
    matm$Cells <- factor(matm$Cells, levels = sort(unique(matm$Cells), decreasing = T))
  }
  matm$colonalID <- factor(matm$colonalID, levels = rev(rownames(data)))

  # Setting color. Lightgrey was Assigned to 0.
  color <- c('lightgrey',RColorBrewer::brewer.pal(n = 9, name = "YlOrRd"))
  names(color) <- c(as.character(0:8),'>8')
  #heatmap
  mytheme <- theme(plot.title = element_blank(), axis.title = element_blank(),
                   axis.ticks = element_blank(),
                   axis.text.y = element_blank(),
                   panel.background = element_blank(),
                   panel.grid = element_blank(),
                   axis.text.x = element_text(angle = 0,vjust = 0.5, hjust = 0.5, size = x.size),
                   legend.title = element_text(family = "ArialMT", size = legend.heatmap.title.size),
                   legend.text = element_text(family = "ArialMT", size = legend.heatmap.size),
                   legend.background = element_rect(fill = "transparent"),
                   plot.margin = unit(c(0,0,0,0), "cm"))
  p <- ggplot(matm, aes(celltype, colonalID)) +
    geom_tile(aes(fill = Cells), colour = NA) +
    scale_fill_manual(values = rev(color)) +
    #scale_x_discrete(position = "top") +
    scale_x_discrete(position = "top",expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    mytheme
  if (ncol(data) > 2) {
    heatleg <- cowplot::get_legend(p)
    hp <- p + theme(legend.position = 'none')
    # Build a legend "bar", row
    legdata$shareTCR <- factor(legdata$shareTCR, levels = unique(legdata$shareTCR))
    legdata$y <- as.numeric(legdata$y)
    leg <- ggplot(legdata, aes(y = y, x = 0)) +
      geom_bar(stat = 'identity', aes(fill = shareTCR)) +
      #geom_point(aes(color = shareTCR), shape = 15, size = 5, show.legend = T) +
      scale_x_discrete(expand = c(1, 0)) + scale_y_discrete(expand = c(0, 0)) +
      scale_fill_manual(values = scPalette2(length(unique(legdata$shareTCR)))) +
      theme_classic() +
      theme(axis.title = element_blank(), axis.line = element_blank(),
            axis.text = element_blank(), axis.ticks = element_blank(),
            legend.title = element_text(family = "ArialMT", size = legend.anno.title.size),
            legend.text = element_text(family = "ArialMT", size = legend.anno.size),
            plot.margin = unit(c(0,0,0,0), "cm"))
    if (!is.null(legend.ncol)) {
      leg <- leg + guides(fill=guide_legend(ncol=legend.ncol))
    }
    annoleg <- cowplot::get_legend(leg)
    anno <- leg + theme(legend.position = 'none')
    #merge legend
    legends <- gridExtra::arrangeGrob(annoleg,heatleg,ncol = 1)
    #merge anno, heatmap and legend
    p <- anno + hp + patchwork::wrap_elements(legends) + patchwork::plot_layout(widths = widths)
  }
  datalist[['plot']] <- p
  if (return.data) {
    return(datalist)
  }else{
    return(p)
  }
}
#' plot TCR heatmap from pre-processed data
#' 
#' @export
plot_TCR_heatmap <- function(data,
                             cell.order = NULL,
                             widths = c(1,8,6), x.size = 15,
                             legend.heatmap.size = 10, 
                             legend.heatmap.title.size = 15,
                             add.legend.anno = T,
                             only.legend.anno = F,
                             legend.ncol = NULL,
                             legend.anno.size = 10,
                             legend.anno.title.size = 15) {
  # generate legend data
  rowprefix <- gsub('\\..*$','', rownames(data))
  shareTCR <- unique(rowprefix)
  y <- as.vector(table(rowprefix)[shareTCR])
  legdata <- data.frame(shareTCR = shareTCR,
                           y = y,
                           stringsAsFactors = F)
  # transform to data frame, matm, for ploting purpose.
  data[data>8] <- 9
  matm <- data.frame(colonalID = rep(rownames(data),ncol(data)),
                     T_celltype = rep(colnames(data),each = nrow(data)),
                     Cells = unlist(as.data.frame(data)))
  matm$celltype <- matm$T_celltype
  if (is.null(cell.order)) {
    cell.order <- unique(matm$celltype)
  }
  matm$celltype <- factor(matm$celltype, levels = cell.order)
  matm$Cells <- as.character(matm$Cells)
  if (any(matm$Cells > 8)) {
    matm$Cells <- gsub('9','>8',matm$Cells)
    matm$Cells <- factor(matm$Cells, levels = rev(c(as.character(0:8),'>8')))
  }else{
    matm$Cells <- factor(matm$Cells, levels = sort(unique(matm$Cells), decreasing = T))
  }
  matm$colonalID <- factor(matm$colonalID, levels = rev(rownames(data)))
  
  # Setting color. Lightgrey was Assigned to 0.
  color <- c('lightgrey',RColorBrewer::brewer.pal(n = 9, name = "YlOrRd"))
  names(color) <- c(as.character(0:8),'>8')
  #heatmap
  mytheme <- theme(plot.title = element_blank(), axis.title = element_blank(),
                   axis.ticks = element_blank(),
                   axis.text.y = element_blank(),
                   panel.background = element_blank(),
                   panel.grid = element_blank(),
                   axis.text.x = element_text(angle = 0,vjust = 0.5, hjust = 0.5, size = x.size),
                   legend.title = element_text(family = "ArialMT", size = legend.heatmap.title.size),
                   legend.text = element_text(family = "ArialMT", size = legend.heatmap.size),
                   legend.background = element_rect(fill = "transparent"),
                   plot.margin = unit(c(0,0,0,0), "cm"))
  p <- ggplot(matm, aes(celltype, colonalID)) +
    geom_tile(aes(fill = Cells), colour = NA) +
    scale_fill_manual(values = rev(color)) +
    #scale_x_discrete(position = "top") +
    scale_x_discrete(position = "top",expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    mytheme
  if (ncol(data) > 2) {
    heatleg <- cowplot::get_legend(p)
    hp <- p + theme(legend.position = 'none')
    # Build a legend "bar", row
    legdata$shareTCR <- factor(legdata$shareTCR, levels = unique(legdata$shareTCR))
    legdata$y <- as.numeric(legdata$y)
    leg <- ggplot(legdata, aes(y = y, x = 0)) +
      geom_bar(stat = 'identity', aes(fill = shareTCR)) +
      #geom_point(aes(color = shareTCR), shape = 15, size = 5, show.legend = T) +
      scale_x_discrete(expand = c(1, 0)) + scale_y_discrete(expand = c(0, 0)) +
      scale_fill_manual(values = scPalette2(length(unique(legdata$shareTCR)))) +
      theme_classic() +
      theme(axis.title = element_blank(), axis.line = element_blank(),
            axis.text = element_blank(), axis.ticks = element_blank(),
            legend.title = element_text(family = "ArialMT", size = legend.anno.title.size),
            legend.text = element_text(family = "ArialMT", size = legend.anno.size),
            plot.margin = unit(c(0,0,0,0), "cm"))
    if (!is.null(legend.ncol)) {
      leg <- leg + guides(fill=guide_legend(ncol=legend.ncol))
    }
    annoleg <- cowplot::get_legend(leg)
    anno <- leg + theme(legend.position = 'none')
    if (only.legend.anno) {
      return(annoleg)
    }
    if (add.legend.anno) {
      #merge legend
      legends <- gridExtra::arrangeGrob(annoleg,heatleg,ncol = 1)
      #merge anno, heatmap and legend
      p <- anno + hp + patchwork::wrap_elements(legends) + patchwork::plot_layout(widths = widths)
    }else{
      p <- anno + hp + heatleg + patchwork::plot_layout(widths = widths)
    }
  }
  return(p)
}
