#' Adjust layer of point to cover outliers.
#' 
#' @return 
#' Ordered point names.
#' @export
cover_outlier <- function(object, annotation, 
                          seurat_clusters = 'seurat_clusters', 
                          reduction = 'umap') {
  df <- cbind(object@reductions[[reduction]]@cell.embeddings[,1:2],
              object@meta.data[,c(annotation, seurat_clusters)])
  colnames(df)[3:4] <- c('annotation','seurat_clusters')
  suerat_cluster_max_anno <- df %>% group_by(seurat_clusters) %>% 
    summarise(max_anno = names(which.max(table(annotation))))
  if (!is.factor(df$annotation)) {
    df$annotation <- factor(df$annotation)
  }
  # add points with high density.
  cell.order <- c()
  for ( i in levels(df$annotation)) {
    cellname <- rownames(df)[df$annotation == suerat_cluster_max_anno$max_anno & 
                               df$seurat_clusters == suerat_cluster_max_anno$seurat_clusters]
    cell.order <- c(cell.order, cellname)
  }
  # add outlier in the end.
  cell.order <- c(cell.order, rownames(df)[!rownames(df) %in% cell.order])
  cell.order <- cell.order[!duplicated(cell.order)]
  return(cell.order)
}
