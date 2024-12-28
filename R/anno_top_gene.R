################## Annotation
#' Annotate raw clusters by max expression, and return an object with annotation Idents
#' @examples
#' markerdf <- data.frame(celltypes = c('T','T','NK','NK','Mono', 'Mono','DC','DC','B','B','Platelet'), 
#' markers = c('CD3D','CD3E','NCAM1','NKG7','CD14','FCGR3A','CST3','CD1C','CD79A','MS4A1','PPBP'))
#' colnames(markerdf)
#' #[1] "celltypes" "markers"
#' pbmc_anno <- anno_top_gene(pbmc, markerdf, group.by = 'seurat_clusters')
#' DimPlot(pbmc_anno, label = T)
#' 
#' @export
anno_top_gene <- function(object, markerdf, group.by = 'seurat_cluster'
) {
  df <- GetAssayData(object)
  df <- df[unique(markerdf$markers),]
  df <- Matrix::t(df)
  ident <- as.vector(object@meta.data[,group.by])
  dflist <- list()
  for (i in unique(ident)) {
    avedf <- apply(df[ident == i,], 2, mean)
    dflist[[paste0('cluster_',i)]] <- avedf
  }
  data <- do.call(rbind, dflist)
  celldflist <- list()
  for (cellt in unique(markerdf$celltypes)) {
    genes <- unique(markerdf$markers[markerdf$celltypes == cellt])
    celldflist[[cellt]] <- apply(data[,genes, drop = F],1,mean)
  }
  celldf <- do.call(cbind, celldflist)
  scaledata <- scale(celldf)
  annotation <- sapply(rownames(scaledata), function(x) colnames(scaledata)[which.max(scaledata[x,])])
  annodf <- data.frame(seurat_clusters = gsub('cluster_','',names(annotation)),
                       Annotation = as.vector(annotation),
                       stringsAsFactors = F)
  colnames(annodf)[1] <- c(group.by)
  meta <- object@meta.data
  meta$Annotation <- NULL
  meta$Cellname <- rownames(meta)
  meta <- left_join(meta, annodf)
  rownames(meta) <- meta$Cellname
  meta$Annotation <- factor(meta$Annotation, levels = unique(markerdf$celltypes))
  Idents(object) <- meta$Annotation
  return(object)
}
#' Annotate raw clusters by max expression, and return an object with annotation Idents
#' @examples
#' markerdf <- data.frame(celltypes = c('T','T','NK','NK','Mono', 'Mono','DC','DC','B','B','Platelet'), 
#' markers = c('CD3D','CD3E','NCAM1','NKG7','CD14','FCGR3A','CST3','CD1C','CD79A','MS4A1','PPBP'))
#' colnames(markerdf)
#' #[1] "celltypes" "markers"
#' pbmc_anno <- anno_top_gene(pbmc, markerdf, group.by = 'seurat_clusters')
#' DimPlot(pbmc_anno, label = T)
#' 
#' @export
anno_top_gene_CODEX <- function(object, markerlist, markerdf = NULL, 
                                group.by = 'seurat_clusters', cut.off = 0.95,
                                assign.noise = F
) {
  df <- GetAssayData(object)
  df <- df[unique(markerlist),]
  df <- Matrix::t(df)
  ident <- as.vector(object@meta.data[,group.by])
  dflist <- list()
  for (i in as.vector(unique(ident))) {
    avedf <- apply(df[ident == i,,drop = F], 2, mean)
    dflist[[paste0('cluster_',i)]] <- avedf
  }
  data <- do.call(rbind, dflist)
  scaledata <- scale(data)
  gene.cut.off <- apply(scaledata, 2, function(x) quantile(x,cut.off))
  annotation <- c()
  annolist <- list()
  for(x in rownames(scaledata)) {
    top_genes <- colnames(scaledata)[scaledata[x,] > gene.cut.off]
    if (length(top_genes) == 0) {
      anno <- 'Noise'
    }else{
      anno <- paste0('Cell_', paste(top_genes, collapse = '_'))
    }
    annotation <- c(annotation,anno)
    annolist[[x]] <- top_genes
  }
  names(annotation) <- rownames(scaledata)
  
  if (!is.null(markerdf)) {
    for (cellt in unique(markerdf$celltypes)) {
      genes <- unique(markerdf$markers[markerdf$celltypes == cellt])
      for (i in 1:length(annolist)) {
        if (sum(genes %in% annolist[[i]]) == length(genes)) annotation[i] <- cellt 
      }
    }
  }
  if (assign.noise) {
    pos <- which(as.vector(annotation) == 'Noise')
    for (i in pos) {
      ave <- c()
      for (cellt in unique(markerdf$celltypes)) {
        genes <- unique(markerdf$markers[markerdf$celltypes == cellt])
        ave <- c(ave, mean(scaledata[i, genes]))
      }
      maxtype <- markerdf$celltypes[which.max(ave)[1]]
      annotation[i] <- maxtype
    }
  }
  annodf <- data.frame(seurat_clusters = gsub('cluster_','',names(annotation)),
                       Annotation = as.vector(annotation),
                       stringsAsFactors = F)
  colnames(annodf)[1] <- c(group.by)
  
  meta <- object@meta.data
  meta$Annotation <- NULL
  meta$Cellname <- rownames(meta)
  meta <- left_join(meta, annodf)
  rownames(meta) <- meta$Cellname
  #meta$Annotation <- factor(meta$Annotation, levels = unique(markerlist))
  Idents(object) <- meta$Annotation
  return(object)
}

anno_top_gene_CODEX_abs <- function(object, markerdf, add.meta = F, gene.cut.off = NULL
) {
  df <- GetAssayData(object)
  df <- df[unique(markerdf$gene),]
  df <- Matrix::t(df)
  
  scaledata <- scale(df)
  
  if (is.null(gene.cut.off)) gene.cut.off <- rep(1, nrow(markerdf))
  
  annotation <- rep('Noise', nrow(df))
  for (i in 1:ncol(df)) {
    pos <- scaledata[,i] > gene.cut.off[i]
    annotation[pos] <- markerdf$celltype[i]
  }
  Idents(object) <- annotation
  
  if (add.meta) {
    for (i in 1:ncol(df)) {
      pos <- scaledata[, i] > gene.cut.off[i]
      object@meta.data[, colnames(df)[i]] <- FALSE
      object@meta.data[pos, colnames(df)[i]] <- TRUE
    }
  }
  
  return(object)
}
