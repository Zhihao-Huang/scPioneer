#' DESeq2 method in Seurat FindMarkers
#' 
#' @export
DEG_Seurat_DESeq2 <- function(object, ident.1, ident.2, 
                              group.by = NULL, 
                              min.pct = 0, logfc.threshold = 0, 
                              verbose = T, p.adjust.method = NULL, ...) {
  DEG <- FindMarkers(object, ident.1 = ident.1, ident.2 = ident.2, 
                     group.by = group.by,  test.use = 'DESeq2', 
                     min.pct = min.pct, logfc.threshold = logfc.threshold, 
                     verbose = verbose, ...)
  DEG <- DEG[order(DEG$avg_log2FC, decreasing = T), ]
  DEG$gene <- rownames(DEG)
  DEG$cluster <- paste0(ident.1, "_vs_", ident.2)
  colnames(DEG)[5] <- "p_val_adj"
  if (!is.null(p.adjust.method)) DEG$p_val_adj <- p.adjust(DEG$p_val, method = p.adjust.method)
  return(DEG)
}