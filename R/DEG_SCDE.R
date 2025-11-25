#' Pairwise DEG using SCDE
#' 
#' @param batchname Colname including covariate in meta.data.
#' @param groupname Colname including two groups for DEG analysis, e.g. SNP and WT labels.
#' @param min.lib.size Minimum number of genes detected in a cell. Cells with fewer genes will be removed (default: 1e3)
#' @param min.reads Minimum number of reads per gene. Genes with fewer reads will be removed (default: 1)
#' @param min.detected Minimum number of cells a gene must be seen in. Genes not seen in a sufficient number of cells will be removed (default: 1)
#' @param length.out number of points (resolution) of the expression magnitude grid (default: 400). Note: larger numbers will linearly increase memory/CPU demands.
#' @param n.randomizations number of bootstrap randomizations to be performed.
#' @param ... Parameter to GetAssayData.
#' 
#' @references https://hms-dbmi.github.io/scde/index.html; https://hms-dbmi.github.io/scw/differential-expression.html
#' @examples
#' pbmc <- readRDS('../data/pbmc_v4.rds')
#' pbmc <- subset(pbmc, Annotation %in% c('B','NK'))
#' pbmc$Annotation <- factor(pbmc$Annotation, levels = unique(pbmc$Annotation))
#' # includes random effects
#' pbmc$batch <- c(rep('sample1',100), rep('sample2',100),rep('sample3',299))
#' DEG <- DEG_SCDE(pbmc, groupname = 'Annotation', ident.1 = 'B', ident.2 = 'NK', batchname  = 'batch', verbose = 0)
#' 
#' @export
DEG_SCDE <- function(object, groupname, ident.1, ident.2 = NULL, batchname = NULL,
                     min.pct = 0.2, logfc.threshold = 0,
                     adjust.method = c('bonferroni','holm', 'hochberg', 'hommel', 'BH', 'BY', 'fdr', 'none'),
                     length.out = 400,
                     n.randomizations  =  100,
                     n.cores = 1, verbose = 1, ...) {
  adjust.method <- match.arg(arg = NULL, choices = adjust.method)
  fc_results <- filter_gene(object = object, group.by = groupname,
                            ident.1 = ident.1, ident.2 = ident.2, min.pct = min.pct, 
                            logfc.threshold = logfc.threshold, slot = 'data', 
                            return.fcresult = T)
  fc_results$gene <- rownames(fc_results)
  # filter genes
  object <- object[rownames(fc_results),]
  # Select groups
  object <- subset(object, cells = colnames(object)[object@meta.data[,groupname] %in% c(ident.1, ident.2)])
  # load raw counts
  cd <- as.matrix(GetAssayData(object,slot = 'counts',...))
  # factor determining cell types
  sg <- object@meta.data[,groupname]
  if (is.null(levels(sg))) {
    sg <- factor(sg, levels = unique(sg))
  }
  # the group factor should be named accordingly
  names(sg) <- colnames(cd)  
  mode(cd) <- 'integer'
  sg <- sg[colnames(cd)]
  message(Sys.time(), ' Calculating SCDE error models using n.cores: ', n.cores,'...')
  o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = n.cores, 
                             min.size.entries = nrow(cd),
                             threshold.segmentation = TRUE, save.crossfit.plots = FALSE,
                             save.model.plots = FALSE, verbose = verbose)
  # filter out cells that don't show positive correlation with
  # the expected expression magnitudes (very poor fits)
  valid.cells <- o.ifm$corr.a > 0
  print(table(valid.cells))
  o.ifm <- o.ifm[valid.cells, ]
  message(Sys.time(), ' Estimating gene expression prior... ')
  o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = length.out, show.plot = FALSE)
  # run differential expression tests on all genes.
  # groups
  groups <- sg[rownames(o.ifm)]
  message(Sys.time(), ' Testing for differential expression...')
  if ( is.null(batchname) | is.numeric(pbmc@meta.data[,batchname]) ) {
    if ( is.numeric(pbmc@meta.data[,batchname]) ) message('Numeric batch is not allowed. Skip batch correction.')
    ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups = groups,
                                        n.randomizations = n.randomizations,
                                        n.cores = n.cores, verbose = verbose)
  }else{
    batch <- object@meta.data[rownames(o.ifm),batchname]
    if (!is.numeric(batch) & is.null(levels(batch))) batch <- factor(batch, levels = unique(batch))
    names(batch) <- rownames(o.ifm)
    ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups = groups,
                                        n.randomizations = n.randomizations,
                                        n.cores = n.cores, verbose = verbose,
                                        batch = batch)
    ediff <- ediff$batch.adjusted
  }
  ediff <- ediff[order(ediff$ce,decreasing = T),]
  ediff$p_val <- 2*pnorm(abs(ediff$Z),lower.tail=F) # 2-tailed p-value
  #ediff$p_val_adj <- 2*pnorm(abs(ediff$cZ),lower.tail=F) # Adjusted to control for FDR
  ediff$p_val_adj <- p.adjust(ediff$p_val, method = adjust.method)
  ediff$gene <- rownames(ediff)
  DEG <- left_join(ediff[, c('p_val','p_val_adj', 'gene')], fc_results, by = 'gene')
  DEG <- DEG[,c('p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj','gene')]
  DEG <- DEG %>% mutate(cluster = ifelse(avg_log2FC < 0, 2, 1)) %>% arrange(cluster,p_val_adj) %>% mutate(cluster = ifelse(cluster == 1, ident.1, ident.2))
  message(Sys.time(), ' Testing for differential expression was done.')
  return(DEG)
}


