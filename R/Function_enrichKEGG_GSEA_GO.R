# clusterProfiler in R 4.3.1 
# source activate easylook

library(clusterProfiler)
library(dplyr)
# over-representation analysis (ORA)

#' GO over-representation analysis
#' 
#' @export
clusterGO <- function (markerlist, 
                       gene_colname = 'gene', cluster_colname = 'cluster',
                       log2FC_colname = 'avg_log2FC',
                       use.top.marker = 100, filter.gene.pattern = NULL, 
                       ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, 
                       geneType = "SYMBOL", OrgDb = "org.Hs.eg.db", n.cores = 1, 
                       ...) 
{
  markerlist <- markerlist[,c(gene_colname, cluster_colname, log2FC_colname)]
  colnames(markerlist) <- c('gene','cluster','avg_log2FC')
  if (!is.null(filter.gene.pattern)) {
    topn <- markerlist %>% filter(!grepl(filter.gene.pattern, 
                                         markerlist$gene)) %>% group_by(cluster) %>% top_n(n = use.top.marker, 
                                                                                           wt = avg_log2FC)
  }
  else {
    topn <- markerlist %>% group_by(cluster) %>% top_n(n = use.top.marker, 
                                                       wt = avg_log2FC)
  }
  GOlist <- parallel::mclapply(unique(topn$cluster), function(t) {
    markers <- topn$gene[topn$cluster == t]
    eg <- clusterProfiler::bitr(markers, fromType = "SYMBOL", 
                                toType = c("ENTREZID", "ENSEMBL",'SYMBOL'), OrgDb = OrgDb)
    #genelist <- unique(eg[, geneType])
    genelist <- markers
    cluster0_markers_go <- clusterProfiler::enrichGO(genelist, 
                                                     OrgDb = OrgDb, ont = ont, pAdjustMethod = pAdjustMethod, 
                                                     pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, 
                                                     keyType = geneType,...)
    df <- as.data.frame(cluster0_markers_go)
    if (nrow(df) == 0) {
      message(paste0("Warning: NO significant GO term found for cluster: ", 
                     t))
      df <- as.data.frame(matrix(nrow = nrow(df), ncol = ncol(df) + 
                                   1, dimnames = list(NULL, c(colnames(df), "celltype"))))
    }
    else {
      df$celltype <- t
    }
    for (i in 1:nrow(df)) {
      geneset <- strsplit(df$geneID[i],'/')[[1]]
      df$geneID[i] <- paste(eg[eg$ENTREZID %in% geneset, 'SYMBOL'], collapse = '/')
    }
    df
  }, mc.cores = n.cores)
  names(GOlist) <- unique(topn$cluster)
  GOmat <- do.call(rbind, GOlist)
  return(GOmat)
}

#' KEGG pathway over-representation analysis
#' 
#' @export
clusterKEGG <- function (markerlist, 
                       gene_colname = 'gene', 
                       log2FC_colname = 'avg_log2FC',
                       log2FC_cutoff = 2,
                       filter.gene.pattern = NULL, 
                       geneType = "SYMBOL", OrgDb = "org.Hs.eg.db", 
                       organism = "hsa", keyType = "kegg", 
                       pvalueCutoff = 0.05,pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                       qvalueCutoff = 0.2, use_internal_data = FALSE,
                       ...) 
{
  markerlist <- markerlist[,c(gene_colname, log2FC_colname)]
  colnames(markerlist) <- c('gene','avg_log2FC')
  if (!is.null(filter.gene.pattern)) {
    topn <- markerlist %>% filter(!grepl(filter.gene.pattern, 
                                         markerlist$gene)) %>% filter(abs(avg_log2FC) > log2FC_cutoff)
  }else{
    topn <- markerlist %>% filter(abs(avg_log2FC) > log2FC_cutoff)
  }
  if (geneType != 'ENTREZID') {
    eg <- clusterProfiler::bitr(markerlist$gene, fromType = geneType, 
                                toType = c("ENTREZID","SYMBOL",'ENSEMBL'), OrgDb = OrgDb)
    genelist <- unique(eg[, 'ENTREZID'])
  }else{
    genelist <- unique(topn$gene)
  }
  
  markers_KEGG <- clusterProfiler::enrichKEGG(genelist, pAdjustMethod = pAdjustMethod, 
                                              pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, 
                                              keyType = keyType,organism = organism,
                                              minGSSize = minGSSize, maxGSSize = maxGSSize,
                                              use_internal_data = use_internal_data, ...)
  df <- as.data.frame(markers_KEGG)
  if (nrow(df) == 0) {
    message("Warning: NO significant KEGG term found")
    return(df)
  }
  for (i in 1:nrow(df)) {
    geneset <- strsplit(df$geneID[i],'/')[[1]]
    df$geneID[i] <- paste(eg[eg$ENTREZID %in% geneset, 'SYMBOL'], collapse = '/')
  }
  return(df)
}

#'  Gene set enrichment analysis
#' 
#' @export
clusterGseKEGG <- function (markerlist, 
                         gene_colname = 'gene', 
                         log2FC_colname = 'avg_log2FC',
                         filter.gene.pattern = NULL, 
                         geneType = "SYMBOL", OrgDb = "org.Hs.eg.db", 
                         organism = 'hsa',
                         minGSSize = 10,
                         maxGSSize = 500,
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH",
                         verbose = TRUE,
                         seed = FALSE,
                         return.gsea_obj = FALSE,
                         by = "fgsea",
                         ...) 
{
  markerlist <- markerlist[,c(gene_colname, log2FC_colname)]
  colnames(markerlist) <- c('gene','avg_log2FC')
  if (!is.null(filter.gene.pattern)) {
    markerlist <- markerlist %>% filter(!grepl(filter.gene.pattern, 
                                               markerlist$gene))
  }
  
  if (geneType != 'ENTREZID') {
    eg <- clusterProfiler::bitr(markerlist$gene, fromType = geneType, 
                                toType = c("ENTREZID",geneType), OrgDb = OrgDb)
    eg$gene <- eg[,geneType]
    markerlist <- left_join(markerlist, eg)
    genelist <- markerlist$avg_log2FC
    names(genelist) <- markerlist$ENTREZID
    genelist <- sort(genelist, decreasing = T)
  }else{
    genelist <- markerlist$avg_log2FC
    names(genelist) <- markerlist$gene
    genelist <- sort(genelist, decreasing = T)
  }

  kegg_gsea <- gseKEGG(genelist,
                       organism = organism,
                       minGSSize = minGSSize,
                       maxGSSize = maxGSSize,
                       pvalueCutoff = pvalueCutoff,
                       pAdjustMethod = pAdjustMethod,
                       verbose = verbose,
                       seed = seed,
                       by = by,
                          ...)
  if (return.gsea_obj) return(kegg_gsea)
  df <- as.data.frame(kegg_gsea)
  if (nrow(df) == 0) {
    message("Warning: NO significant term found")
    return(df)
  }
  for (i in 1:nrow(df)) {
    if (is.na(df$core_enrichment[i])) next
    geneset <- strsplit(df$core_enrichment[i],'/')[[1]]
    df$core_enrichment[i] <- paste(eg[eg$ENTREZID %in% geneset, 'SYMBOL'], collapse = '/')
  }
  return(df)
}

#'  Gene set enrichment analysis
#' 
#' @export
clusterGESA <- function (markerlist, 
                         gene_colname = 'gene', 
                         log2FC_colname = 'avg_log2FC',
                         filter.gene.pattern = NULL, 
                         geneType = "SYMBOL", OrgDb = "org.Hs.eg.db", 
                         exponent = 1,
                         minGSSize = 10,
                         maxGSSize = 500,
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH",
                         verbose = TRUE,
                         seed = FALSE,
                         return.gsea_obj = FALSE,
                         by = "fgsea",
                         ...) 
{
  markerlist <- markerlist[,c(gene_colname, log2FC_colname)]
  colnames(markerlist) <- c('gene','avg_log2FC')
  if (!is.null(filter.gene.pattern)) {
    markerlist <- markerlist %>% filter(!grepl(filter.gene.pattern, 
                                         markerlist$gene))
  }
  
  if (geneType != 'ENTREZID') {
    eg <- clusterProfiler::bitr(markerlist$gene, fromType = geneType, 
                                toType = c("ENTREZID",geneType), OrgDb = OrgDb)
    eg$gene <- eg[,geneType]
    markerlist <- left_join(markerlist, eg)
    genelist <- markerlist$avg_log2FC
    names(genelist) <- markerlist$ENTREZID
    genelist <- sort(genelist, decreasing = T)
  }else{
    genelist <- markerlist$avg_log2FC
    names(genelist) <- markerlist$gene
    genelist <- sort(genelist, decreasing = T)
  }
  
  markers_gsea <- DOSE::gseDO(genelist,
                              exponent = exponent,
                              minGSSize = minGSSize,
                              maxGSSize = maxGSSize,
                              pvalueCutoff = pvalueCutoff,
                              pAdjustMethod = pAdjustMethod,
                              verbose = verbose,
                              seed = seed,
                              by = by,
                              ...)
  if (return.gsea_obj) return(markers_gsea)
  df <- as.data.frame(markers_gsea)
  if (nrow(df) == 0) {
    message("Warning: NO significant term found")
    return(df)
  }
  for (i in 1:nrow(df)) {
    if (is.na(df$core_enrichment[i])) next
    geneset <- strsplit(df$core_enrichment[i],'/')[[1]]
    df$core_enrichment[i] <- paste(eg[eg$ENTREZID %in% geneset, 'SYMBOL'], collapse = '/')
  }
  return(df)
}
