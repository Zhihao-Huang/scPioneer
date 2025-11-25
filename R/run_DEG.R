#' DEG methods for benchmarking single-cell RNA-seq data.
#' 
#' 
#' @param object A Seurat object.
#' @param groupname Column name of meta.data including 2 groups to performing differential expression.
#' @param ident.1 One cluster in groupname.
#' @param ident.2 Another cluster in groupname.
#' @param covariates Column names in meta.data for confounding effect correction.
#' @param min.pct Minimum fraction of min.pct cells in either of the two populations in groupname.
#' @param logfc.threshold Minimum log2 Fold-change value
#' @param test.method test method to use.
#' @param pseudobulk Whether to perform pseudobulk anlaysis by Seurat::AggregateExpression.
#' @param pseudobulk_covar covariates for pseudobulk generation.
#' @param pseudo_count pseudo count for log normalization in some method (DESeq2).
#' @param adjust.method Adjusted method for P-value.
#' @param n.cores Number of cores for parallel processes.
#' @param workdir Output directory for temporary result.
#' @param ... Arguments passed to other methods and to specific DE methods. 
#' @examples 
#' pbmc <- readRDS('../data/pbmc_v4.rds')
#' pbmc <- subset(pbmc, Annotation %in% c('B','NK'))
#' pbmc$Annotation <- factor(pbmc$Annotation, levels = unique(pbmc$Annotation))
#' # includes random effects
#' pbmc$batch <- c(rep('sample1',100), rep('sample2',100),rep('sample3',299))
#' resultlist <- run_DEG(pbmc, 'Annotation', 'B', 'NK', covariates = 'batch', test.method = c('Seurat_t','glmGamPoi'), min.pct = 0.3, logfc.threshold = 1)
#' resultlistp <- run_DEG(pbmc, 'Annotation', 'B', 'NK', pseudobulk = T, test.method = c('Seurat_wilcox','glmGamPoi'), pseudobulk_covar = 'batch',min.pct = 0.3, logfc.threshold = 1)
#' @export
run_DEG <- function(object, groupname, ident.1, ident.2 = NULL, 
                    covariates = NULL,
                    min.pct = 0, logfc.threshold = 0,
                    min.cell = 50,
                    test.method = c('Seurat_wilcox', 'Seurat_t', 'Seurat_negbinom', 
                                    'Seurat_poisson', 'Seurat_DESeq2',
                                    't_test', 'poisson_Wald', 'negbinom_Wald',
                                    'MAST', 'scvi',
                                    'limma_trend','limma_voom',
                                    'glmGamPoi','DESeq2_LRT', 'DESeq2_Wald', ### time consuming , not recommended
                                    'edgeR_LRT', 'edgeR_QLF',
                                    'Monocle3_quasipoisson', 'Monocle3_negbinomial',
                                    'Monocle3_poisson', 
                                    'Monocle3_gaussian', 'Monocle3_zipoisson',
                                    'Monocle3_zinegbinomial',
                                    'Monocle3_mixed.negbinomial',
                                    'Monocle3_binomial', ## return 0 genes, not recommended.
                                    'SCDE','miloDE',
                                    'DAseq_STG','DAseq_wilcox', 'DAseq_t', 'DAseq_negbinom', 
                                    'DAseq_poisson', 'DAseq_MAST', 'DAseq_DESeq2',
                                    'Meld_rank', 'Meld_t-test', 'Meld_wald', 'Meld_lrt',
                                    'Meld_Seurat_wilcox', 'Meld_Seurat_t', 'Meld_Seurat_negbinom', 
                                    'Meld_Seurat_poisson', 'Meld_Seurat_MAST', 'Meld_Seurat_DESeq2'),
                    pseudo_count = 1,
                    adjust.method = c('bonferroni','holm', 'hochberg', 'hommel', 'BH', 'BY', 'fdr', 'none'),
                    n.cores = 1, workdir = NULL, 
                    python2use_Meld = '/media/london_B/huangzh/software/miniconda3/envs/FR-Perturb/bin/python',
                    python2use_scvi = '/lustre/home/kwxiong/Huangzhihao/software/miniconda3/envs/scvi-env/bin/python',
                    meld_script_path = '/lustre/home/kwxiong/Huangzhihao/project/benchmark_perturb/shell/Perturb_benchmark-main/script/Meld_all.py', 
                    meta_script_path = '/lustre/home/kwxiong/Huangzhihao/project/benchmark_perturb/shell/Perturb_benchmark-main/script/Meld_meta.py', 
                    deg_script_path = '/lustre/home/kwxiong/Huangzhihao/project/benchmark_perturb/shell/Perturb_benchmark-main/script/Meld_deg.py',
                    classes_script_path = '/lustre/home/kwxiong/Huangzhihao/project/benchmark_perturb/shell/Perturb_benchmark-main/script/Meld_classes.py',
                    scvi_script_path = '/lustre/home/kwxiong/Huangzhihao/project/benchmark_perturb/shell/Perturb_benchmark-main/script/scvi.py',
                    ...) {
  start.time <- Sys.time()
  adjust.method <- match.arg(arg = NULL, choices = adjust.method)
  # checking groups
  if (any(table(object@meta.data[,groupname]) < min.cell)) {
    message(paste0('Error: one of group has cells fewer than ', min.cell, '. Increase min.cell if you want to continue DE analysis.'))
    return(data.frame())
  }
  
  ## add pseudo_count to avoid error when performing DESeq2
  if (grepl('DESeq2', test.method)) {
    object[["pseudocount"]] <- CreateAssayObject(counts = as.matrix(object[["RNA"]]@counts)+pseudo_count)
    DefaultAssay(object) <- 'pseudocount'
  }
  
  DEGlist <- list()
  Seurat_method_pos <- test.method %in% c('Seurat_wilcox', 'Seurat_t', 'Seurat_negbinom', 
                                          'Seurat_poisson')
  if (any(Seurat_method_pos)) {
    test.use <- test.method[Seurat_method_pos]
    for (i in test.use) {
      message(paste0('##################  Perform ',i,'   ###################'))
      DEG <- do.call(paste0('DEG_', i), list(object, ident.1 = ident.1, ident.2 = ident.2, 
                        group.by = groupname, 
                        min.pct = min.pct, logfc.threshold = logfc.threshold,
                        latent.vars = NULL, p.adjust.method = adjust.method,
                        verbose = T, ...)
                     )
      DEGlist[[paste0('Seurat_',i)]] <- DEG
    }
  }
  if ('Seurat_DESeq2' %in% test.method) {
    message('##################        Perform Seurat_DESeq2         ##################')
    DEG <- FindMarkers(object, ident.1 = ident.1, ident.2 = ident.2,
                       group.by = groupname,  test.use = 'DESeq2',
                       min.pct = min.pct, logfc.threshold = logfc.threshold,
                       verbose = T)
    DEG <- DEG[order(DEG$avg_log2FC, decreasing = T), ]
    DEG$gene <- rownames(DEG)
    DEG$cluster <- paste0(ident.1, "_vs_", ident.2)
    colnames(DEG)[5] <- "p_val_adj"
    if (!is.null(adjust.method)) DEG$p_val_adj <- p.adjust(DEG$p_val, method = adjust.method)
    DEGlist[['Seurat_DESeq2']] <- DEG
  }
  if ('scvi' %in% test.method) {
    message('##################        Perform scvi         ##################')
    DEG <- DEG_scvi(object, ident.1 = ident.1, ident.2 = ident.2, 
                  groupname = groupname, 
                  covariates = covariates, 
                  min.pct = min.pct, logfc.threshold = logfc.threshold,
                  python2use = python2use_scvi,
                  scvi_script_path = scvi_script_path,
                  n.cores = n.cores, ...)
    DEGlist[['scvi']] <- DEG
  }
  if ('t_test' %in% test.method) {
    message('##################    Perform lm (t-test)      ##################')
    DEG <- DEG_t_test(object, ident.1 = ident.1, ident.2 = ident.2, 
                  groupname = groupname, 
                  covariates = covariates, 
                  min.pct = min.pct, logfc.threshold = logfc.threshold,
                  p.adjust.method = adjust.method,
                  n.cores = n.cores, ...)
    DEGlist[['t_test']] <- DEG
  }
  if ('poisson_Wald' %in% test.method) {
    message('######## Perform Wald test based on poisson model ############')
    DEG <- DEG_poisson_Wald(object, ident.1 = ident.1, ident.2 = ident.2, 
                           groupname = groupname, 
                           min.pct = min.pct, logfc.threshold = logfc.threshold,
                           covariates = covariates,
                           p.adjust.method = adjust.method,
                           n.cores = n.cores, ...)
    DEGlist[['poisson_Wald']] <- DEG
  }
  if ('negbinom_Wald' %in% test.method) {
    message('##### Perform Wald test based on negative binomial model ######')
    DEG <- DEG_negbinom_Wald(object, ident.1 = ident.1, ident.2 = ident.2, 
                            groupname = groupname, 
                            min.pct = min.pct, logfc.threshold = logfc.threshold,
                            covariates = covariates,
                            p.adjust.method = adjust.method,
                            n.cores = n.cores, ...)
    DEGlist[['negbinom_Wald']] <- DEG
  }
  if ('MAST' %in% test.method) {
    message('##################        Perform MAST         ##################')
    DEG <- DEG_MAST(object,group.colname = groupname, ident.1 = ident.1, ident.2 = ident.2, 
                    covariates = covariates, pseudo_count = pseudo_count, 
                    min.cell = 0, parallel = T, n.cores = n.cores, 
                    min.pct = min.pct, logfc.threshold = logfc.threshold,
                    p.adjust.method = adjust.method, ...)
    DEGlist[['MAST']] <- DEG
  }
  if ('limma_trend' %in% test.method) {
    message('##################     Perform limma_trend     ##################')
    DEG <- DEG_limma_trend(object, group.colname = groupname, ident.1 = ident.1, ident.2 = ident.2, 
                           batch.colname = covariates, prior.count = pseudo_count, 
                           min.pct = min.pct, logfc.threshold = logfc.threshold,
                           adjust.method = adjust.method, ...)
    DEGlist[['limma_trend']] <- DEG
  }
  if ('limma_voom' %in% test.method) {
    message('##################     Perform limma_voom     ##################')
    DEG <- DEG_limma_voom(object, group.colname = groupname,ident.1 = ident.1, ident.2 = ident.2, 
                          min.pct = min.pct, logfc.threshold = logfc.threshold,
                          batch.colname = covariates, adjust.method = adjust.method, ...)
    DEGlist[['limma_voom']] <- DEG
  }
  if ('DESeq2_LRT' %in% test.method) {
    message('##################    Perform DESeq2_LRT     ##################')
    DEG <- DEG_DESeq2(object, groupname = groupname,  ident.1 = ident.1, ident.2 = ident.2,
                      adjust.method = adjust.method,
                      fitType = 'glmGamPoi', useT = T, alpha = 0.9999,
                      cooksCutoff = Inf,
                      min.pct = min.pct, logfc.threshold = logfc.threshold,
                      independentFiltering = F,
                      test.method = 'LRT', batchname  = covariates,
                      n.cores = n.cores, ...)
    DEGlist[['DESeq2_LRT']] <- DEG
  }
  if ('DESeq2_Wald' %in% test.method) {
    message('##################    Perform DESeq2_Wald    ##################')
    DEG <- DEG_DESeq2(object, groupname = groupname, ident.1 = ident.1, ident.2 = ident.2,
                      adjust.method = adjust.method, alpha = 0.9999, 
                      cooksCutoff = Inf,
                      independentFiltering = F,
                      min.pct = min.pct, logfc.threshold = logfc.threshold,
                      test.method = 'Wald', parallel = T, n.cores = n.cores, useT = F,
                      fitType = 'local', batchname  = covariates, ...)
    DEGlist[['DESeq2_Wald']] <- DEG
  }
  if ('edgeR_LRT' %in% test.method) {
    message('##################     Perform edgeR_LRT     ##################')
    DEG <- DEG_edgeR(object, groupname = groupname,
                     ident.1 = ident.1, ident.2 = ident.2, adjust.method = adjust.method,
                     min.pct = min.pct, logfc.threshold = logfc.threshold,
                     test = 'LRT', batchname  = covariates, ...)
    DEGlist[['edgeR_LRT']] <- DEG
  }
  if ('edgeR_QLF' %in% test.method) {
    message('##################    Perform edgeR_QLF      ##################')
    DEG <- DEG_edgeR(object, groupname = groupname,
                     ident.1 = ident.1, ident.2 = ident.2, adjust.method = adjust.method,
                     min.pct = min.pct, logfc.threshold = logfc.threshold,
                     test = 'QLF', batchname  = covariates, ...)
    DEGlist[['edgeR_QLF']] <- DEG
  }
  if (grepl('^Monocle3',test.method)) {
    Monocle3.method <- gsub('^Monocle3_','',test.method[grepl('^Monocle3',test.method)])
    for (i in Monocle3.method) {
      i <- gsub('\\.','-',i)
      message(paste0('##################  Perform Monocle3 ',i,'  #################'))
      DEG <- DEG_Monocle3(object = object, groupname = groupname, ident.1 = ident.1, ident.2 = ident.2,
                          batchname = covariates, adjust.method = adjust.method,
                          min.pct = min.pct, logfc.threshold = logfc.threshold,
                          expression_family = i,
                          n.cores = n.cores, ...)
      DEGlist[[paste0('Monocle3_', i)]] <- DEG
    }
  }
  if ('SCDE' %in% test.method) {
    message('##################       Perform SCDE        ##################')
    DEG <- DEG_SCDE(object, groupname = groupname, #batchname  = covariates, 
                    ident.1 = ident.1, ident.2 = ident.2,
                    length.out = 400, adjust.method = adjust.method,
                    min.pct = min.pct, logfc.threshold = logfc.threshold,
                    n.randomizations  =  100,
                    n.cores = n.cores, verbose = 0, ...)
    DEGlist[['SCDE']] <- DEG
  }
  if ('miloDE' %in% test.method) {
    message('##################       Perform miloDE        ##################')
    DEG <- DEG_miloDE(object, 
                      groupname = groupname, ident.1 = ident.1, ident.2 = ident.2,
                      covariates = covariates,
                      min.pct = min.pct, logfc.threshold = logfc.threshold,
                      min.cell = min.cell, 
                      n.cores = n.cores, ...)
    DEGlist[['miloDE']] <- DEG
  }
  if (grepl('DAseq',test.method)) {
    DAseq.method <- gsub('DAseq_','',test.method[grepl('DAseq',test.method)])
    for (i in DAseq.method) {
      message(paste0('##################  Perform DAseq ',i,'  #################'))
      DEG <- DEG_DAseq(object, groupname = groupname, p.adjust.method = adjust.method,
                       min.pct = min.pct, logfc.threshold = logfc.threshold,
                       ident.1 = ident.1, ident.2 = ident.2, 
                       min.cell = min.cell, reduction = 'pca',
                       python2use = python2use_Meld,
                       ...)
      DEGlist[[paste0('DAseq_',i)]] <- DEG
    }
  }
  if (grepl('Meld',test.method)) {
    Meld.method <- gsub('Meld_','',test.method[grepl('Meld',test.method)])
    for (i in Meld.method) {
      message(paste0('##################  Perform Meld ',i,'  #################'))
      test.Seurat = NULL
      if(grepl('Seurat',i)) {
        test.Seurat <- gsub('^.*Seurat_','',i)
        test.meld <- 'rank'
      }else{
        test.Seurat <- NULL
        test.meld <- i
      }
      tmpdir <- ifelse(is.null(workdir), tempdir(), workdir)
      DEG <- DEG_Meld(object, groupname, ident.1, ident.2, p.adjust.method = adjust.method,
                      covariates = covariates, test.Seurat = test.Seurat,
                      min.pct = min.pct, logfc.threshold = logfc.threshold,
                      min.cell = min.cell, 
                      test.meld = test.meld, tmpdir = tmpdir,
                      python2use = python2use_Meld, 
                      meld_script_path = meld_script_path,
                      meta_script_path = meta_script_path,
                      deg_script_path = deg_script_path,
                      classes_script_path = classes_script_path,
                      ...
      )
      DEGlist[[paste0('Meld_',i)]] <- DEG
    }
  }
  
  message('########   Processing ',paste(c('Group',ident.1, test.method), collapse = ', '),' done.  #########')
  end.time <- Sys.time()
  time.use <- difftime(end.time, start.time, units = "secs")
  message(paste0('Time usage: ', round(time.use, 4), ' seconds.'))
  DEGlist[['time.mins']] <- time.use
  return(DEGlist)
}
