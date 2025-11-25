#' Pairwise DEG using limma-trend
#' # citation: https://ucdavis-bioinformatics-training.github.io/2019-single-cell-RNA-sequencing-Workshop-UCD_UCSF/scrnaseq_analysis/scRNA_Workshop-PART6.html
#' # Refer to https://github.com/eonurk/RNA-seq-differential-analyses-guideline
#' 
#' @param batch.colname Colname including covariate in meta.data.
#' @param group.colname Colname including two groups for DEG analysis, e.g. SNP and WT labels.
#' @param min.cell Genes with number of expressed cells smaller than min.cell will be removed.
#' @param max.cpm.cut.off Filter low-expressed genes.
#' @param calcNormFactors Calculate normalization factors from counts data.
#' @param ... Parameter to GetAssayData.
#' @examples
#' pbmc <- readRDS('../data/pbmc_v4.rds')
#' DEG <- DEG_limma_trend(pbmc, 'Annotation', ident.1 = 'Memory CD4 T', ident.2 = 'Naive CD4 T')
#' # includes random effects
#' pbmc$batch <- c(rep('sample1',1000), rep('sample2',1000),rep('sample3',638))
#' DEG <- DEG_limma_trend(pbmc, 'Annotation', ident.1 = 'Memory CD4 T', ident.2 = 'Naive CD4 T', batch.colname = c('nCount_RNA','batch'))
#' @export
DEG_limma_trend <- function(object, group.colname, ident.1, ident.2,  
                            batch.colname = NULL, 
                            min.pct = 0.2, logfc.threshold = 0,
                            calcNormFactors = T, 
                            prior.count = 1,
                            duplicateCorrelation.name = NULL, 
                            trend = T,
                            adjust.method = c('bonferroni','holm', 'hochberg', 'hommel', 'BH', 'BY', 'fdr', 'none'),
                            sort.by = "P",...) {
  adjust.method <- match.arg(arg = NULL, choices = adjust.method)
  # filter genes
  fc_results <- filter_gene(object = object, group.by = group.colname,
                            ident.1 = ident.1, ident.2 = ident.2, min.pct = min.pct, 
                            logfc.threshold = logfc.threshold, slot = 'data', 
                            return.fcresult = T)
  fc_results$gene <- rownames(fc_results)
  object <- object[rownames(fc_results),]
  # Perform CPM normalization or get normalized data directly from object.
  if (calcNormFactors) {
    counts <- as.matrix(GetAssayData(object, slot = 'counts', assay = 'RNA'))
    d <- DGEList(counts)
    d <- calcNormFactors(d)
    expr <- edgeR::cpm(d, log=TRUE, prior.count = prior.count)
  }else{
    expr <- as.matrix(GetAssayData(object,...))
  }
  # create disign matrix
  object$group.colname <- object@meta.data[,group.colname]
  object$group.colname <- as.vector(object$group.colname)
  object$group.colname <- check_replace_name(object$group.colname,
                                             '\\+| ',  '.', fixed = F)
  ident.1.adj <- check_replace_name(ident.1, '\\+| ',  '.', fixed = F)
  ident.2.adj <- check_replace_name(ident.2, '\\+| ',  '.', fixed = F)
  if (is.null(batch.colname)) {
    mm <- model.matrix(~0 + group.colname, data = object@meta.data)
  }else{
    mdf <- object@meta.data[, c('group.colname'), drop = F]
    for (i in batch.colname) {
      if (is.numeric(object@meta.data[,i])) {
        mdf[, paste0('covariate_', i)] <- object@meta.data[,i]
      }else{
        mdf[, paste0('covariate_', i)] <- check_replace_name(as.vector(object@meta.data[,i]), '\\+| ',  '.', fixed = F)
      }
    }
    mm <- model.matrix(~0 + ., data = mdf)
  }
  # Estimate the intra-block correlation given a block structure for the arrays or samples.
  if (!is.null(duplicateCorrelation.name)) {
    object$dupC.name <- object@meta.data[,duplicateCorrelation.name]
    object$dupC.name <- as.vector(object$dupC.name)
    object$dupC.name <- check_replace_name(object$dupC.name,
                                           '\\+| ',  '.', fixed = F)
    corfit <- duplicateCorrelation(expr, mm, block=object$dupC.name)
    fit <- lmFit(expr, mm, block=object$dupC.name, correlation=corfit$consensus)
  }else{
    fit <- lmFit(expr, mm)
  }
  x <- paste0('group.colname', ident.1.adj,' - ', 'group.colname', ident.2.adj)
  contr <- makeContrasts(contrasts = x, levels = colnames(coef(fit)))
  DEG <- contrasts.fit(fit, contrasts = contr)
  # trend is TRUE in this function.
  DEG <- eBayes(DEG, trend=trend)
  top.table <- topTable(DEG, sort.by = sort.by, n = Inf, adjust.method = adjust.method)
  top.table$gene <- rownames(top.table)
  colnames(top.table)[4] <- 'p_val'
  colnames(top.table)[5] <- 'p_val_adj'
  DEG <- left_join(top.table[, c('p_val','p_val_adj', 'gene')], fc_results, by = 'gene')
  DEG <- DEG[,c('p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj','gene')]
  DEG <- DEG %>% mutate(cluster = ifelse(avg_log2FC < 0, 2, 1)) %>% arrange(cluster,p_val_adj) %>% mutate(cluster = ifelse(cluster == 1, ident.1, ident.2))
  return(DEG)
}