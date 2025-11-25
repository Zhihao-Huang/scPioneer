#' Pairwise DEG using DEseq2 adjusted for single-cell data.
#' Refer to https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#recommendations-for-single-cell-analysis. https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/05b_wald_test_results.html
#' 
#' @param test Likelihood ratio test (LRT) or Wald test.
#' @param fitType Parameter in function DESeq(). Either "parametric", "local", "mean", or "glmGamPoi" for the type of fitting of dispersions to the mean intensity. See estimateDispersions for description.
#' @param useT Parameter in function DESeq(). Logical, passed to nbinomWaldTest, default is FALSE, where Wald statistics are assumed to follow a standard Normal.
#' @param minReplicatesForReplace Parameter in function DESeq(). The minimum number of replicates required in order to use replaceOutliers on a sample. If there are samples with so many replicates, the model will be refit after these replacing outliers, flagged by Cook's distance. Set to Inf in order to never replace outliers. It set to Inf for fitType="glmGamPoi".
#' @param minmu Lower bound on the estimated count for fitting gene-wise dispersion and for usage with nbinomWaldTest and nbinomLRT. If fitType="glmGamPoi", then 1e-6 will be used (as this fitType is optimized for single cell data, where a lower minmu is recommended), otherwise the default value as evaluated on bulk datasets is 0.5.
#' @param ... Parameter to GetAssayData.
#' 
#' @examples
#' pbmc <- readRDS('../data/pbmc_v4.rds')
#' pbmc <- subset(pbmc, Annotation %in% c('B','NK'))
#' pbmc$Annotation <- factor(pbmc$Annotation, levels = unique(pbmc$Annotation))
#' # includes random effects
#' pbmc$batch <- c(rep('sample1',100), rep('sample2',100),rep('sample3',299))
#' DEG_LRT <- DEG_DESeq2(pbmc, groupname = 'Annotation', useT = T, test.method = 'LRT', batchname  = 'batch')
#' DEG_Wald <- DEG_DESeq2(pbmc, groupname = 'Annotation', test.method = 'Wald', parallel = T, fitType = 'parametric', batchname  = 'batch')
#' 
#' @export
DEG_DESeq2 <- function(object, groupname, ident.1, ident.2,
                       batchname = NULL, assay = 'RNA',
                       min.pct = 0.2, logfc.threshold = 0,
                       filtergene.rowSums = 0,
                       test.method = c("LRT","Wald"),
                       fitType = c("glmGamPoi"),
                       useT = F,
                       minReplicatesForReplace = Inf,
                       calucateSumFactors = F,
                       minmu = if (fitType == "glmGamPoi") 1e-06 else 0.5,
                       adjust.method = c('bonferroni','holm', 'hochberg', 'hommel', 'BH', 'BY', 'fdr', 'none'),
                       parallel = F,
                       n.cores = 4,
                       alpha = 0.9999,
                       cooksCutoff = Inf,
                       independentFiltering = F,
                       ...) {
  test.method <- match.arg(arg = NULL, choices = test.method)
  adjust.method <- match.arg(arg = NULL, choices = adjust.method)
  fitType <- match.arg(arg = NULL, choices = fitType)
  # filter genes
  fc_results <- filter_gene(object = object, group.by = groupname,
                            ident.1 = ident.1, ident.2 = ident.2, min.pct = min.pct, 
                            logfc.threshold = logfc.threshold, slot = 'data', assay = assay,
                            return.fcresult = T)
  fc_results$gene <- rownames(fc_results)
  object <- object[rownames(fc_results),]
  # design matrix
  object$group.colname <- object@meta.data[,groupname]
  cellID <- colnames(object)[object@meta.data[,groupname] %in% c(ident.1, ident.2)]
  object <- subset(object, cells = cellID)
  object$group.colname <- factor(object$group.colname, levels = c(ident.1, ident.2))
  level <- levels(object$group.colname)
  object$group.colname <- as.vector(object$group.colname)
  object$group.colname <- check_replace_name(object$group.colname, '\\+| ',  '.', fixed = F)
  object$group.colname <- factor(object$group.colname, levels = check_replace_name(level, '\\+| ',  '.', fixed = F))

  
  if (is.null(batchname)) {
    mdf <- object@meta.data[, c('group.colname'), drop = F]
    reduced <- as.formula( "~1")
    form <- ~group.colname
  }else{
    mdf <- object@meta.data[, c('group.colname'), drop = F]
    for (i in batchname) {
      if (is.numeric(object@meta.data[,i])) {
        mdf[, paste0('covariate_', i)] <- object@meta.data[,i]
      }else{
        mdf[, paste0('covariate_', i)] <- check_replace_name(as.vector(object@meta.data[,i]), '\\+| ',  '.', fixed = F)
        mdf[, paste0('covariate_', i)] <- factor(mdf[, paste0('covariate_', i)], levels = unique(mdf[, paste0('covariate_', i)]))  
      }
    }
    reduced <- as.formula(paste( "~", base::paste(paste0('covariate_',batchname), sep = ' + ')))
    form = as.formula(paste( "~", base::paste('group.colname',paste0('covariate_',batchname), sep = ' + ')))
  }
  counts <- GetAssayData(object, layer = 'counts',...)
  # Create DESeq object
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = mdf,
                                design= ~ group.colname)
  # to avoid error in calculating log geometric means
  dds <- DESeq(dds, sfType = "poscounts")
  # Renew design
  design(dds) <- form
  # Reorder levels
  dds$group.colname <- relevel(dds$group.colname, ref = ident.2)
  # filter genes by rowSums
  dds <- dds[rowSums(counts(dds)) > filtergene.rowSums,]
  if (calucateSumFactors) dds <- scran::computeSumFactors(dds, BPPARAM = MulticoreParam(n.cores))
  
  if (test.method == 'Wald') {
    message('Performing DESeq-Wald test...')
    dds <- DESeq(dds, test = test.method, useT = useT, fitType = fitType, 
                 minReplicatesForReplace = minReplicatesForReplace,
                 parallel = parallel, minmu = minmu,
                 BPPARAM = MulticoreParam(n.cores)
    )
  }else{
    message('Performing DESeq-LRT test...')
    dds <- DESeq(dds, test = test.method, useT = useT, fitType = fitType, 
                 reduced = reduced, 
                 minReplicatesForReplace = minReplicatesForReplace,
                 parallel = parallel, minmu = minmu,
                 BPPARAM = MulticoreParam(n.cores)
    )
  }
  
  print(resultsNames(dds)) # lists the coefficients
  resultgname <- resultsNames(dds)[grepl('^group.colname',resultsNames(dds))]
  print(resultgname)
  resultgname <- resultgname[1]
  res <- results(dds, name = resultgname, 
                 pAdjustMethod = adjust.method,
                 alpha = alpha,
                 cooksCutoff = cooksCutoff,
                 independentFiltering = independentFiltering,
                 parallel = parallel,
                 BPPARAM = MulticoreParam(n.cores))
  res$gene <- rownames(res)
  colnames(res)[5] <- 'p_val'
  colnames(res)[6] <- 'p_val_adj'
  res <- as.data.frame(res)
  DEG <- left_join(res[, c('p_val','p_val_adj', 'gene')], fc_results, by = 'gene')
  DEG <- DEG %>% mutate(cluster = ifelse(avg_log2FC < 0, 2, 1)) %>% arrange(cluster,p_val_adj) %>% mutate(cluster = ifelse(cluster == 1, ident.1, ident.2))
  return(DEG)
}