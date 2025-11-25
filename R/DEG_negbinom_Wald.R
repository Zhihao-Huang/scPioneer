#' Fitting gene expression using Negative Binomial Generalized linear model. 
#' 
#' @param object A Seurat object.
#' @param groupname Column name of meta.data including 2 groups to performing differential expression.
#' @param ident.1 One cluster in groupname.
#' @param ident.2 Another cluster in groupname.
#' @param min.cell The minimum number of cells in each group.
#' @param covariates Column names in meta.data for confounding effect correction.
#' 
#' @examples
#' pbmc <- scPioneer::pbmc
#' pbmc$batch <- c(rep('sample1',100), rep('sample2',100),rep('sample3',299))
#' coefdf <- DEG_negbinom_Wald(pbmc, 'Annotation', 'B','NK', n.cores = 8)
#' coefdf <- DEG_negbinom_Wald(pbmc, 'Annotation', 'B','NK', covariates = c('batch','percent.mt'), n.cores = 8)
#' 
#' @references  https://stat.ethz.ch/R-manual/R-devel/library/MASS/html/glm.nb.html
#' @export
DEG_negbinom_Wald <- function(object, groupname, ident.1, ident.2, p.adjust.method = 'fdr',
                             min.pct = 0.2, logfc.threshold = 0,
                             covariates = NULL, n.cores = 4 ) {
  object <- subset(object, cells = colnames(object)[object@meta.data[,groupname] %in% c(ident.1, ident.2)])
  fc_results <- filter_gene(object = object, group.by = groupname,
                            ident.1 = ident.1, ident.2 = ident.2, min.pct = min.pct, 
                            logfc.threshold = logfc.threshold, slot = 'data', 
                            return.fcresult = T)
  # filter genes
  object <- object[rownames(fc_results),]
  # Fit linear model for each gene.
  message(paste0('Fitting a Negative Binomial Generalized linear model ...'))
  coeflist <- mclapply(rownames(object), function(x) {
    # Create covariate data frame
    sample_df <- data.frame(values = GetAssayData(object, slot = 'counts')[x,])
    sample_df$groups <- factor(object@meta.data[, groupname], levels = c(ident.1, ident.2))
    if (!is.null(covariates)) sample_df <- cbind(sample_df, object@meta.data[, covariates])
    sample_df <- as.data.frame(sample_df)
    rownames(sample_df) <- colnames(object)
    colnames(sample_df) <- c('values','groups', covariates)
    # Create formula
    formula_str <- paste0("values ~ ", paste(c('groups', covariates), collapse = " + "))
    fmla <- as.formula(formula_str)
    # change reference level - groups1 in output is groups with 2nd level, NK, as base
    fit.glm <- MASS::glm.nb(fmla, 
                            data = sample_df,
                            contrasts = list(groups = contr.treatment(2, base = 2)))
    coef_result <- coef(summary(fit.glm))
    return(coef_result[2,])
  }, mc.cores = n.cores)
  # merge results of all genes 
  coefdf <- do.call(rbind, coeflist)
  coefdf <- as.data.frame(coefdf)
  if (dim( coefdf )[1] == 0) return(coefdf)
  rownames(coefdf) <- rownames(object)
  colnames(coefdf)[4] <- 'p_val'
  coefdf$p_val_adj <- p.adjust(coefdf$p_val, method = p.adjust.method)
  # Ordered by p_val_adj
  coefdf$gene <- rownames(coefdf)
  fc_results$gene <- rownames(fc_results)
  DEG <- left_join(coefdf, fc_results, by = 'gene')
  DEG <- DEG[,c('p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj','gene')]
  DEG <- DEG %>% mutate(cluster = ifelse(avg_log2FC < 0, 2, 1)) %>% arrange(cluster,p_val_adj) %>% mutate(cluster = ifelse(cluster == 1, ident.1, ident.2))
  
  return(DEG)
}

