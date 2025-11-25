#' Pairwise DEG using limma-voom.
#' Refer to https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
#' 
#' @param batch.colname Colname including covariate in meta.data.
#' @param group.colname Colname including two groups for DEG analysis, e.g. SNP and WT labels.
#' @param min.cell Genes with number of expressed cells smaller than min.cell will be removed.
#' @param max.cpm.cut.off Filter low-expressed genes.
#' @param ... Parameter to GetAssayData.
#' @examples
#' pbmc <- readRDS('../data/pbmc_v4.rds')
#' DEG <- DEG_limma_voom(pbmc, 'Annotation', ident.1 = 'Memory CD4 T', ident.2 = 'Naive CD4 T')
#' # includes random effects
#' pbmc$batch <- c(rep('sample1',1000), rep('sample2',1000),rep('sample3',638))
#' DEG <- DEG_limma_voom(pbmc, 'Annotation', ident.1 = 'Memory CD4 T', ident.2 = 'Naive CD4 T', batch.colname = 'batch')
#' @export
DEG_limma_voom <- function(object, group.colname, ident.1, ident.2,  
                           min.pct = 0.2, logfc.threshold = 0,
                           batch.colname = NULL, 
                           adjust.method = c('bonferroni','holm', 'hochberg', 'hommel', 'BH', 'BY', 'fdr', 'none'),
                           duplicateCorrelation.name = NULL, 
                           plot.voom = F,
                           ...) {
  adjust.method <- match.arg(arg = NULL, choices = adjust.method)
  # filter genes
  fc_results <- filter_gene(object = object, group.by = group.colname,
                            ident.1 = ident.1, ident.2 = ident.2, min.pct = min.pct, 
                            logfc.threshold = logfc.threshold, slot = 'data', 
                            return.fcresult = T)
  fc_results$gene <- rownames(fc_results)
  object <- object[rownames(fc_results),]
  
  counts <- as.matrix(GetAssayData(object, slot = 'counts', ...))
  DEGobject <- DGEList(counts)
  DEGobject <- calcNormFactors(DEGobject)
  
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
  if (!is.null(duplicateCorrelation.name)) {
    object$dupC.name <- object@meta.data[,duplicateCorrelation.name]
    object$dupC.name <- as.vector(object$dupC.name)
    object$dupC.name <- check_replace_name(object$dupC.name,
                                           '\\+| ',  '.', fixed = F)
    corfit <- duplicateCorrelation(expr, mm, block=object$dupC.name)
    y <- voom(DEGobject, mm, plot = plot.voom)
    fit <- lmFit(y, mm, block=object$dupC.name, correlation=corfit$consensus)
  }else{
    y <- voom(DEGobject, mm, plot = plot.voom)
    fit <- lmFit(y, mm) 
  }
  x <- paste0('group.colname', ident.1.adj,' - ', 'group.colname', ident.2.adj)
  contr <- makeContrasts(contrasts = x, levels = colnames(coef(fit)))
  DEG <- contrasts.fit(fit, contrasts = contr)
  DEG <- eBayes(DEG)
  top.table <- topTable(DEG, sort.by = "P", n = Inf, adjust.method = adjust.method)
  top.table$gene <- rownames(top.table)
  colnames(top.table)[4] <- 'p_val'
  colnames(top.table)[5] <- 'p_val_adj'
  DEG <- left_join(top.table[, c('p_val','p_val_adj', 'gene')], fc_results, by = 'gene')
  DEG <- DEG[,c('p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj','gene')]
  DEG <- DEG %>% mutate(cluster = ifelse(avg_log2FC < 0, 2, 1)) %>% arrange(cluster,p_val_adj) %>% mutate(cluster = ifelse(cluster == 1, ident.1, ident.2))
  return(DEG)
}