#' Differential expression using MAST
#' 
#' @param object A Seurat object.
#' @param groupname Column name of meta.data including 2 groups to performing differential expression.
#' @param ident.1 One cluster in groupname.
#' @param ident.2 Another cluster in groupname.
#' @param min.cell The minimum number of cells in each group.
#' @param covariates Column names in meta.data for confounding effect correction.
#' @param random_effect.vars Covariates accounted for random effect.
#' @param method Argument to zlm. Character vector, either 'glm', 'glmer' or 'bayesglm'.
#' @param ebayes Argument to zlm. If TRUE, regularize variance using empirical bayes method.
#' @param n.cores How many cores to use.
#' 
#' @references Seurat::MASTDETest; Andrew McDavid, Greg Finak and Masanao Yajima (2017). 
#' MAST: Model-based Analysis of Single Cell Transcriptomics. R package version 1.26.0. 
#' https://github.com/RGLab/MAST/. https://github.com/satijalab/seurat/issues/3712
#' 
#' @examples
#' pbmc <- scPioneer::pbmc
#' pbmc$batch <- c(rep('sample1',1000), rep('sample2',1000),rep('sample3',638))
#' subs4 <- subset(pbmc, Annotation %in% c('CD14+ Mono','FCGR3A+ Mono'))
#' subs4 <- subs4[rowSums(GetAssayData(subs4, slot = 'counts')) > 100, ]
#' MASTdf <- DEG_MAST(subs4, 'Annotation', ident.1 = 'CD14+ Mono', ident.2 = 'FCGR3A+ Mono', covariates = 'percent.mt')
#' MASTdf <- DEG_MAST(subs4, 'Annotation', ident.1 = 'CD14+ Mono', ident.2 = 'FCGR3A+ Mono', covariates = 'percent.mt', random_effect.vars = 'batch', method = 'glmer', ebayes = F)
#' 
#' @export
DEG_MAST <- function(object, group.colname, ident.1, ident.2,  
                     covariates = NULL, random_effect.vars = NULL,
                     min.pct = 0.2, logfc.threshold = 0,
                     method	= 'bayesglm', ebayes = T,
                     pseudo_count = 1,
                     parallel = T, n.cores = 4, 
                     p.adjust.method = 'fdr',...
) {
  op <- options("mc.cores"=n.cores)
  on.exit(options(op), add = T)
  
  object <- subset(object, cells = colnames(object)[object@meta.data[,group.colname] %in% c(ident.1, ident.2)])
  fc_results <- filter_gene(object = object, group.by = group.colname,
                            ident.1 = ident.1, ident.2 = ident.2, min.pct = min.pct, 
                            logfc.threshold = logfc.threshold, slot = 'data', 
                            return.fcresult = T)
  fc_results$gene <- rownames(fc_results)
  # filter genes
  object <- object[rownames(fc_results),]
  # Cell metadata
  coldata <- object@meta.data[,c(group.colname, covariates, random_effect.vars), drop = F]
  # For auto-ordering DEG result, change names of elements in groups.
  coldata$Group <- 'ident1'
  coldata$Group[object@meta.data[,group.colname] == ident.2] <- 'ident2'
  coldata <- coldata[,c('Group', covariates, random_effect.vars), drop = F]
  # Gene metadata
  fData <- data.frame(gene=rownames(object))
  # Create SingleCellAssay object
  sca <- suppressMessages(MAST::FromMatrix(exprsArray=log2(as.matrix(GetAssayData(object, slot = 'counts', assay = 'RNA')) + pseudo_count),
                                           cData=coldata, 
                                           fData=fData))
  SummarizedExperiment::colData(sca)$Group <- relevel(factor(SummarizedExperiment::colData(sca)$Group), ref="ident2")
  
  message(paste0('Fitting model using ',method,' ...'))
  
  if (!is.null(random_effect.vars)){
    formula_str <- paste0(" ~ ", paste(c('Group', covariates), collapse = " + "),
                          paste(" + (1 |", random_effect.vars, ")", collapse=""))
    message(formula_str)
    fmla <- as.formula(formula_str)
    # fit model with random effect. Recommended parameters: method = 'glmer', ebayes = F,strictConvergence = FALSE
    zlmCond <- MAST::zlm(formula = fmla, sca = sca, method = method,
                         ebayes = ebayes, strictConvergence = F, ...)
  }else{
    formula_str <- paste0(" ~ ", paste(c('Group', covariates), collapse = " + "))
    message(formula_str)
    fmla <- as.formula(formula_str)
    zlmCond <- MAST::zlm(formula = fmla, sca = sca, method= method, 
                         ebayes = ebayes, parallel = parallel, ...)
  }
  
  message('Summarizing ...')
  contrast_name <- paste0('Group', 'ident1')
  summaryCond <- suppressMessages(MAST::summary(zlmCond, doLRT = contrast_name))
  summaryDt <- summaryCond$datatable
  fcHurdle <- merge(summaryDt[summaryDt$contrast==contrast_name
                              & summaryDt$component=='logFC', c(1,7,5,6,8)],
                    summaryDt[summaryDt$contrast==contrast_name
                              & summaryDt$component=='H', c(1,4)],
                    by = 'primerid')
  
  fcHurdle <- stats::na.omit(as.data.frame(fcHurdle))
  colnames(fcHurdle)[1] <- 'gene'
  colnames(fcHurdle)[6] <- 'p_val'
  
  if (nrow(fcHurdle) < 1) return(fcHurdle)
  
  fcHurdle$p_val_adj <- p.adjust(fcHurdle$p_val, p.adjust.method)
  DEG <- left_join(fcHurdle[, c('p_val','p_val_adj', 'gene')], fc_results, by = 'gene')
  DEG <- DEG[,c('p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj','gene')]
  DEG <- DEG %>% mutate(cluster = ifelse(avg_log2FC < 0, 2, 1)) %>% arrange(cluster,p_val_adj) %>% mutate(cluster = ifelse(cluster == 1, ident.1, ident.2))
  
  return(DEG)
}
