#' Pairwise DEG using edegeR. 
#' QLF: Conduct genewise quasi F-tests for a given coefficient or coefficient contrast.
#' LRT: conduct likelihood ratio tests for a given coefficient or coefficient contrast.
#' 
#' @param ... Parameter to GetAssayData.
#' @param test likelihood ratio test (LRT) or quasi-likelihood F-test (QLF)
#' @param min.cells the minimum number of expressed cells for each genes.
#' @examples
#' pbmc <- readRDS('../data/pbmc_v4.rds')
#' pbmc <- subset(pbmc, Annotation %in% c('B','NK'))
#' pbmc$Annotation <- factor(pbmc$Annotation, levels = unique(pbmc$Annotation))
#' # includes other covariates
#' pbmc$batch <- c(rep('sample1',100), rep('sample2',100),rep('sample3',299))
#' DEG <- DEG_edgeR(pbmc, groupname = 'Annotation', ident.1 = 'B', ident.2 = 'NK', batchname  = c('batch','percent.mt'))
#' DEG <- DEG_edgeR(pbmc, groupname = 'Annotation', ident.1 = 'B', ident.2 = 'NK', batchname  = c('batch','percent.mt'), test = 'QLF')
#' @export
DEG_edgeR <- function(object, groupname, ident.1, ident.2, batchname = NULL,
                      min.pct = 0.2, logfc.threshold = 0,
                      test = c("LRT",'QLF'),
                      adjust.method = c('bonferroni','holm', 'hochberg', 'hommel', 'BH', 'BY', 'fdr', 'none'),
                      ...) {
  test <- match.arg(arg = NULL, choices = test)
  adjust.method <- match.arg(arg = NULL, choices = adjust.method)
  fc_results <- filter_gene(object = object, group.by = groupname,
                            ident.1 = ident.1, ident.2 = ident.2, min.pct = min.pct, 
                            logfc.threshold = logfc.threshold, slot = 'data', 
                            return.fcresult = T)
  fc_results$gene <- rownames(fc_results)
  # filter genes
  object <- object[rownames(fc_results),]
  
  # trimming names
  object$group.colname <- as.vector(object@meta.data[, groupname])
  object$group.colname <- check_replace_name(object$group.colname,
                                             '\\+| ',  '.', fixed = F)
  ident.1.adj <- check_replace_name(ident.1, '\\+| ',  '.', fixed = F)
  ident.2.adj <- check_replace_name(ident.2, '\\+| ',  '.', fixed = F)
  object$group.colname <- factor(object$group.colname, levels = c(ident.2.adj,ident.1.adj))
  # create DGEList
  counts_mat <- GetAssayData(object, slot = 'counts',...)
  dge <- DGEList(
    counts = counts_mat, 
    group = object$group.colname
  )
  # design matrix
  mdf <- object@meta.data[, c('group.colname'), drop = F]
  if (!is.null(batchname)) {
    for (i in batchname) {
      if (is.numeric(object@meta.data[,i])) {
        mdf[, paste0('covariate_', i)] <- object@meta.data[,i]
      }else{
        mdf[, paste0('covariate_', i)] <- check_replace_name(as.vector(object@meta.data[,i]), '\\+| ',  '.', fixed = F)
        mdf[, paste0('covariate_', i)] <- factor(mdf[, paste0('covariate_', i)], levels = unique(mdf[, paste0('covariate_', i)]))
      }
    }
  }
  design <- model.matrix(~., data = mdf)

  # fit models
  dge <- calcNormFactors(dge)
  if (test == 'QLF') {
    dge <- estimateDisp(dge, design = design)
    fit <- glmQLFit(dge, design)
    res <- glmQLFTest(fit, coef = 2)
    deg_QLF <- topTags(res,n = nrow(res), adjust.method = adjust.method, sort.by = 'logFC')
    DEG <- deg_QLF$table
    DEG$cluster <- paste0(ident.1.adj,'_vs_',ident.2.adj)
  }else if (test == 'LRT') {
    dge <- estimateGLMCommonDisp(dge,design)
    dge <- estimateGLMTrendedDisp(dge,design, method="auto")
    dge <- estimateGLMTagwiseDisp(dge,design)
    fit <- glmFit(dge, design)
    lrt12 <- glmLRT(fit, coef = 2)
    deg_LRT <- topTags(lrt12,n = nrow(lrt12), adjust.method = adjust.method, sort.by = 'logFC')
    DEG <- deg_LRT$table
    DEG$cluster <- paste0(levels(mdf$group.colname)[2],'_vs_',levels(mdf$group.colname)[1])
  }
  colnames(DEG)[4] <- 'p_val'
  colnames(DEG)[5] <- 'p_val_adj'
  DEG$gene <- rownames(DEG)
  DEG <- left_join(DEG[, c('p_val','p_val_adj', 'gene')], fc_results, by = 'gene')
  DEG <- DEG[,c('p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj','gene')]
  DEG <- DEG %>% mutate(cluster = ifelse(avg_log2FC < 0, 2, 1)) %>% arrange(cluster,p_val_adj) %>% mutate(cluster = ifelse(cluster == 1, ident.1, ident.2))
  
  return(DEG)
}


