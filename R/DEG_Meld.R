#' Pairwise DEG using Meld
#' Refer to https://nbviewer.org/github/KrishnaswamyLab/MELD/blob/main/notebooks/MELD_thresholding.Tcell.ipynb
#' 
#' @param on_disk Whether to save temporary results on disk.
#' @param test.Seurat Employ DE methods built in Seurat.
#' @param test.meld Employ DE methods built in Meld.
#' @param return.classes whether to return clusters of cells.
#' @param ... Parameter to GetAssayData.
#' 
#' @examples
#' pbmc <- readRDS('../data/pbmc_v4.rds')
#' pbmc <- subset(pbmc, Annotation %in% c('B','NK'))
#' pbmc$Annotation <- factor(pbmc$Annotation, levels = unique(pbmc$Annotation))
#' # includes random effects
#' pbmc$batch <- c(rep('sample1',100), rep('sample2',100),rep('sample3',299))
#' DEG <- DEG_Meld(pbmc, 'Annotation', 'B','NK')
#' DEG <- DEG_Meld(pbmc, 'Annotation', 'B','NK', test.Seurat = 'MAST', covariates = c('batch', 'percent.mt'))
#' 
#' @export
DEG_Meld <- function(object, groupname, ident.1, ident.2,
                     covariates = NULL, test.Seurat = NULL,
                     python2use = '/media/london_B/huangzh/software/miniconda3/envs/FR-Perturb/bin/python',
                     meld_script_path = './Meld_all.py',
                     meta_script_path = './Meld_meta.py',
                     deg_script_path = './Meld_deg.py',
                     classes_script_path = './Meld_classes.py',
                     min.cell = 50,
                     min.pct = 0.2, logfc.threshold = 0,
                     on_disk = T, tmpdir = tempdir(),
                     test.meld = c('rank', 't-test','wald','lrt'),
                     p.adjust.method = c('bonferroni','holm', 'hochberg', 'hommel', 'BH', 'BY', 'fdr', 'none'),
                     return.classes = F
) {
  test.meld <- match.arg(arg = NULL, choices = test.meld) 
  p.adjust.method <- match.arg(arg = NULL, choices = p.adjust.method)
  assignInNamespace("is_conda_python", function(x){ return(FALSE) }, ns="reticulate")
  if (!is.null(python2use)) reticulate::use_python(python2use)
  print(reticulate::py_config())
  if(!suppressWarnings(reticulate::py_module_available("meld"))) stop("python module meld does not seem to be installed; - try running 'py_config()'")
  
  reticulate::source_python(meld_script_path)
  reticulate::source_python(meta_script_path)
  reticulate::source_python(deg_script_path)
  
  # filter genes
  fc_results <- filter_gene(object = object, group.by = groupname,
                            ident.1 = ident.1, ident.2 = ident.2, min.pct = min.pct, 
                            logfc.threshold = logfc.threshold, slot = 'data', 
                            return.fcresult = T)
  fc_results$gene <- rownames(fc_results)
  object <- object[rownames(fc_results),]
  #### Input expression matrix
  mat <- GetAssayData(object, slot = 'data')
  print(class(mat))
  #X <- bigdataTOmat(mat, min_num = 20000, transpose = T, as.type = "TsparseMatrix")
  X <- as(Matrix::t(mat), "TsparseMatrix")
  i <- as.integer(X@i)
  j <- as.integer(X@j)
  val <- X@x
  dim <- as.integer(X@Dim)
  #### Input meta data
  group <- as.vector(object@meta.data[,groupname])
  group <- gsub(ident.1,'stimulated', group)
  group <- gsub(ident.2,'unstimulated', group)
  DEG_Meld_py_args<-c(list(i=i, j=j, val=val, dim=dim, cell_idx=colnames(object),
                           feature_idx = rownames(object),
                           group=group, test_method=test.meld))
  if (on_disk) {
    # Memory-no-mapping error occurs when calculating Gaussian Mixture by reticulate-python. 
    # Employ python script for Calculating GaussianMixture instead of using R package reticulate.
    message('Calculating likehood metadata...')
    metadata <- do.call(meta_Meld_py, DEG_Meld_py_args)
    filename <- paste0(ident.1, '_vs_',ident.2)
    message(filename)
    metapath <-paste0(tmpdir, '/meld_meta_', filename, '.csv')
    write.csv(metadata, file = metapath)
    classpath <- paste0(tmpdir, '/meld_class_', filename, '.csv')
    message('Calculating classes using Gaussian Mixture...')
    system(paste0(python2use,' ',classes_script_path, ' --inputpath ', metapath, ' --outputpath ',classpath))
    classes <- read.csv(classpath, row.names = 1)
    classes <- as.vector(classes[,1])
    if (any(!c(0,2) %in% classes)) {
      message(filename, ' has less than 2 groups in classes. Return empty.')
      return(list(classes = classes, de_results = NULL))
    }
    message('Calculating DEG...')
    if (is.null(test.Seurat)) {
      deg_Meld_py_args <- c(list(i=i, j=j, val=val, dim=dim, cell_idx=colnames(object),
                                 feature_idx = rownames(object),
                                 classes=classes, test_method=test.meld))
      deg <- do.call(deg_Meld_py, deg_Meld_py_args)
      Meld_results <- list(classes = classes, de_results = deg)
    }else{
      Meld_results <- list(classes = classes, de_results = NULL)
    }
  }else{
    Meld_results <- do.call(DEG_Meld_py, DEG_Meld_py_args)
    names(Meld_results) <- c('classes','de_results')
    if (Meld_results$de_results == '') Meld_results$de_results <- NULL
  }
  if (return.classes) return(Meld_results[['classes']])
  # Check cell number of depleted and enriched group in classes.
  if (sum(classes == 0) < min.cell | sum(classes == 2) < min.cell) {
    message(paste0('Number of cells is fewer than ',min.cell, '. Return Empty.'))
    return(list(classes = classes, de_results = data.frame()))
  }
  if (!is.null(test.Seurat)) {
    object$classes <- Meld_results[['classes']]
    print(table(object$classes))
    object <- subset(object, classes %in% c(0,2))
    DEG <- FindMarkers(object, ident.1 = 2, ident.2 = 0, 
                       group.by = 'classes', test.use = test.Seurat, 
                       min.pct = 0, logfc.threshold = 0, 
                       latent.vars = covariates,
                       verbose = T)
    DEG$gene <- rownames(DEG)
    DEG$cluster <- paste0(ident.1, "_vs_", ident.2)
    colnames(DEG)[5] <- "p_val_adj"
    DEG$p_val_adj <- p.adjust(DEG$p_val, method = p.adjust.method)
    DEG <- DEG[order(DEG$p_val_adj),]
    Meld_results[['de_results']] <- DEG
  }else{
    DEG <- as.data.frame(Meld_results[['de_results']])
    print(head(DEG))
    DEG$p_val_adj <- p.adjust(DEG$pval, method = p.adjust.method)
    colnames(DEG)[4] <- 'avg_log2FC'
    DEG$cluster <- paste0(ident.1, "_vs_", ident.2)
    DEG <- DEG[order(DEG$p_val_adj),]
    Meld_results[['de_results']] <- DEG
  }
  if (is.null(Meld_results[['de_results']]$DEG)) return(Meld_results)
  DEG <- left_join(Meld_results[['de_results']]$DEG[, c('p_val','p_val_adj', 'gene')], fc_results, by = 'gene')
  DEG <- DEG[,c('p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj','gene')]
  DEG <- DEG %>% mutate(cluster = ifelse(avg_log2FC < 0, 2, 1)) %>% arrange(cluster,p_val_adj) %>% mutate(cluster = ifelse(cluster == 1, ident.1, ident.2))
  Meld_results[['de_results']]$DEG <- DEG
  return(Meld_results)
}
