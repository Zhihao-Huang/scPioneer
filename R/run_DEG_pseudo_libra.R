#' DEG analysis using libra.
#' 
#' @param object A Seurat object.
#' @param groupname Column name of meta.data including 2 groups to performing differential expression.
#' @param ident.1 One cluster in groupname.
#' @param ident.2 Another cluster in groupname.
#' @param covariates Column names in meta.data for confounding effect correction.
#' @param min.pct Minimum fraction of min.pct cells in either of the two populations in groupname.
#' @param logfc.threshold Minimum log2 Fold-change value
#' @param test.method test method to use.
#' @param de_family DEG family to use. Available options are: singlecell, pseudobulk and mixedmodel. Default is pseudobulk.
#' @param de_method DEG methods to use. Only 3 methods are available for psudobulk analysis: 'DESeq2','edgeR', 'limma'.
#' @param de_type Sub-types of DEG, including LRT, QLF, Wald, trend, voom. 
#' @param min_reps Argument to run_de. The minimum number of replicates in a cell type to retain it. Defaults to 2. description
#' @param covariates latent variables for single-cell Seurat/Signac based methods.
#' @param adjust.method Adjusted method for P-value.
#' @param n.cores Number of cores for parallel processes.
#' @param ... Arguments passed to other methods and to specific DE methods. 
#' @examples 
#' pbmc <- readRDS('../data/pbmc_v4.rds')
#' pbmc <- subset(pbmc, Annotation %in% c('B','NK'))
#' pbmc$Annotation <- factor(pbmc$Annotation, levels = unique(pbmc$Annotation))
#' # includes random effects
#' pbmc$batch <- c(rep('sample1',100), rep('sample2',100),rep('sample3',299))
#' DEG <- run_DEG_pseudo_libra(object = pbmc, celltype = 'orig.ident', samplename = 'batch', groupname = 'Annotation', ident.1 = 'B', ident.2 = 'NK')
#' @export
run_DEG_pseudo_libra <- function(object, groupname, ident.1, ident.2 = NULL, 
                                 samplename = NULL, celltype =NULL,
                                 min.pct = 0, 
                                 min.cells = 50,
                                 min_reps = 2,
                                 min_features = 0,
                                 de_family = 'pseudobulk',
                                 de_method = c('DESeq2','edgeR', 'limma'),
                                 de_type = c('LRT', 'QLF', 'Wald', 'trend', 'voom'),
                                 adjust.method = c('fdr','bonferroni','holm', 'hochberg', 'hommel', 'BH', 'BY',  'none'),
                                 n.cores = 1,  ...) {
  de_method <- match.arg(arg = NULL, choices = de_method)
  de_type <- match.arg(arg = NULL, choices = de_type)
  adjust.method <- match.arg(arg = NULL, choices = adjust.method)
  if (is.null(samplename)) samplename <- 'orig.ident'
  # checking groups
  if (any(table(object@meta.data[,groupname]) < min.cells)) {
    message(paste0('Error: one of group has cells fewer than ', min.cells, '. Increase min.cell if you want to continue DE analysis.'))
    return(data.frame())
  }
  # filtering genes
  fcresult <- FoldChange(object, ident.1, 
                         group.by = groupname, assay = 'RNA', slot = 'counts')
  # feature selection (based on percentages)
  alpha.min <- pmax(fcresult$pct.1, fcresult$pct.2)
  names(x = alpha.min) <- rownames(x = fcresult)
  features <- names(x = which(x = alpha.min >= min.pct))
  genemeta <- fcresult[features,]
  message(paste0('Number of genes remained: ', nrow(genemeta)))
  object <- object[rownames(genemeta),]
  object@meta.data[,groupname] <- factor(object@meta.data[,groupname], levels = c(ident.1, ident.2))
  
  # pesudo-bulk DEG
  DE = run_de(input = object, replicate_col = samplename, cell_type_col = celltype, label_col = groupname,
              de_family = de_family, de_method = de_method, de_type = de_type,
              n_threads = n.cores,
              min_cells = min.cells,
              min_reps = min_reps,
              min_features = min_features,...)
  DE <- DE[order(DE$avg_logFC,decreasing = T),]
  DE <- DE %>% group_by(cell_type) %>% mutate(p_val_adj = p.adjust(p_val, method = adjust.method))
  return(DE)
}
