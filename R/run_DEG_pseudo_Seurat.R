#' Check whether there is special symbol in names of cell barcode or annotation to avoid error.
#'
#' @param celltype A vector, can be barcode or annotation.
#'
#' @export
check_replace_name <- function(celltype, pattern = NULL, replacement = NULL, fixed = F) {
  # user pattern
  if (!is.null(pattern) & !is.null(replacement)) {
    if (any(grepl(pattern = pattern, x = celltype, fixed = fixed))) {
      warning(paste0("Feature names cannot have characters (",pattern,"), replacing with ampersand (",replacement,")"),
              call. = FALSE, immediate. = TRUE)
      celltype <- gsub(pattern = pattern, replacement = replacement,
                       x = celltype, fixed = fixed)
      
    }
  }
  return(celltype)
}
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
#' resultlistp <- run_DEG(pbmc, 'Annotation', 'B', 'NK', pseudobulk = T, test.method = c('Seurat_wilcox','glmGamPoi'), pseudobulk_covar = 'batch',min.pct = 0.3, logfc.threshold = 1)
#' @export
run_DEG_pseudo_Seurat <- function(object, groupname, ident.1, ident.2 = NULL, 
                    min.pct = 0, logfc.threshold = 0, 
                    min.cell = 50,
                    test.method = c('Seurat_wilcox', 'Seurat_t', 'Seurat_negbinom', 
                                    'Seurat_poisson', 'Seurat_MAST', 'Seurat_DESeq2',
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
                    pseudobulk_covar = NULL, 
                    pseudo_count = 1,
                    adjust.method = c('bonferroni','holm', 'hochberg', 'hommel', 'BH', 'BY', 'fdr', 'none'),
                    n.cores = 1, workdir = NULL, ...) {
  adjust.method <- match.arg(arg = NULL, choices = adjust.method)
  # checking groups
  if (any(table(object@meta.data[,groupname]) < min.cell)) {
    message(paste0('Error: one of group has cells fewer than ', min.cell, '. Increase min.cell if you want to continue DE analysis.'))
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
  genemeta <- genemeta[abs(genemeta$avg_log2FC) >= logfc.threshold, ]
  message(paste0('Number of genes remained: ', nrow(genemeta)))
  object <- object[rownames(genemeta),]
  object@meta.data[,groupname] <- factor(object@meta.data[,groupname], levels = c(ident.2, ident.1))
  
  # pesudo-bulk
  # pseudobulk the counts based on donor-condition

  object <- AggregateExpression(object, assays = "RNA", return.seurat = T, group.by = c(groupname, pseudobulk_covar))
  Idents(object) <- object@meta.data[, groupname]


  ## add pseudo_count to avoid error when performing DESeq2
  if (grepl('DESeq2', test.method)) {
    object[["pseudocount"]] <- CreateAssayObject(counts = as.matrix(object[["RNA"]]@counts)+pseudo_count)
    DefaultAssay(object) <- 'pseudocount'
  }
  
  DEGlist <- list()
  Seurat_method_pos <- test.method %in% c('Seurat_wilcox', 'Seurat_t', 'Seurat_negbinom', 
                                          'Seurat_poisson', 'Seurat_MAST', 'Seurat_DESeq2')
  if (any(Seurat_method_pos)) {
    test.use <- gsub('Seurat_','',test.method[Seurat_method_pos])
    for (i in test.use) {
      message(paste0('##################  Perform Seurat_',i,'   ###################'))
      DEG <- do.call(paste0('DEG_Seurat_', i), list(object, ident.1 = ident.1, ident.2 = ident.2, 
                                                    group.by = groupname, test.use = i, 
                                                    min.pct = 0, logfc.threshold = 0, 
                                                    latent.vars = NULL,
                                                    verbose = T, ...)
      )
      DEG <- DEG[order(DEG$avg_log2FC, decreasing = T), ]
      DEG$gene <- rownames(DEG)
      DEG$cluster <- paste0(ident.1, "_vs_", ident.2)
      colnames(DEG)[5] <- "p_val_adj"
      DEGlist[[paste0('Seurat_',i)]] <- DEG
    }
  }
  if ('limma_trend' %in% test.method) {
    message('##################     Perform limma_trend     ##################')
    DEG <- DEG_limma_trend(object, group.colname = groupname, ident.1 = ident.1, ident.2 = ident.2, 
                           batch.colname = covariates, prior.count = pseudo_count, 
                           adjust.method = adjust.method, ...)
    colnames(DEG)[4] <- 'p_val'
    colnames(DEG)[5] <- 'p_val_adj'
    #colnames(DEG)[1] <- 'avg_log2FC'
    DEGlist[['limma_trend']] <- DEG
  }
  if ('limma_voom' %in% test.method) {
    message('##################     Perform limma_voom     ##################')
    DEG <- DEG_limma_voom(object, group.colname = groupname,ident.1 = ident.1, ident.2 = ident.2, 
                          batch.colname = covariates, adjust.method = adjust.method, ...)
    colnames(DEG)[4] <- 'p_val'
    colnames(DEG)[5] <- 'p_val_adj'
    # colnames(DEG)[1] <- 'avg_log2FC'
    DEGlist[['limma_voom']] <- DEG
  }
  if ('DESeq2_LRT' %in% test.method) {
    message('##################    Perform DESeq2_LRT     ##################')
    DEG <- DEG_DESeq2(object, groupname = groupname, 
                      adjust.method = adjust.method,
                      fitType = 'glmGamPoi', useT = T, alpha = 0.9999,
                      cooksCutoff = Inf,
                      independentFiltering = F,
                      test.method = 'LRT', batchname  = covariates,
                      n.cores = n.cores, ...)
    DEG$gene <- rownames(DEG)
    #colnames(DEG)[2] <- 'avg_log2FC'
    colnames(DEG)[5] <- 'p_val'
    colnames(DEG)[6] <- 'p_val_adj'
    DEGlist[['DESeq2_LRT']] <- DEG
  }
  if ('DESeq2_Wald' %in% test.method) {
    message('##################    Perform DESeq2_Wald    ##################')
    DEG <- DEG_DESeq2(object, groupname = groupname, 
                      adjust.method = adjust.method, alpha = 0.9999, 
                      cooksCutoff = Inf,
                      independentFiltering = F,
                      test.method = 'Wald', parallel = T, n.cores = n.cores, useT = F,
                      fitType = 'local', batchname  = covariates, ...)
    DEG$gene <- rownames(DEG)
    #colnames(DEG)[2] <- 'avg_log2FC'
    colnames(DEG)[5] <- 'p_val'
    colnames(DEG)[6] <- 'p_val_adj'
    DEGlist[['DESeq2_Wald']] <- DEG
  }
  if ('edgeR_LRT' %in% test.method) {
    message('##################     Perform edgeR_LRT     ##################')
    DEG <- DEG_edgeR(object, groupname = groupname,
                     ident.1 = ident.1, ident.2 = ident.2, adjust.method = adjust.method,
                     test = 'LRT', batchname  = covariates, ...)
    DEG$gene <- rownames(DEG)
    #colnames(DEG)[1] <- 'avg_log2FC'
    colnames(DEG)[4] <- 'p_val'
    colnames(DEG)[5] <- 'p_val_adj'
    DEGlist[['edgeR_LRT']] <- DEG
  }
  if ('edgeR_QLF' %in% test.method) {
    message('##################    Perform edgeR_QLF      ##################')
    DEG <- DEG_edgeR(object, groupname = groupname,
                     ident.1 = ident.1, ident.2 = ident.2, adjust.method = adjust.method,
                     test = 'QLF', batchname  = covariates, ...)
    DEG$gene <- rownames(DEG)
    #colnames(DEG)[1] <- 'avg_log2FC'
    colnames(DEG)[4] <- 'p_val'
    colnames(DEG)[5] <- 'p_val_adj'
    DEG <- DEG[order(DEG$logFC, decreasing = T),]
    DEGlist[['edgeR_QLF']] <- DEG
  }
  if (grepl('^Monocle3',test.method)) {
    Monocle3.method <- gsub('^Monocle3','',test.method[grepl('^Monocle3',test.method)])
    for (i in Monocle3.method) {
      i <- gsub('\\.','-',i)
      message(paste0('##################  Perform Monocle3 ',i,'  #################'))
      DEG <- DEG_Monocle3(object = object, groupname = groupname, 
                          batchname = covariates, adjust.method = adjust.method,
                          expression_family = i,
                          n.cores = n.cores, ...)
      # DEG$cluster <- paste0(ident.1, "_vs_", ident.2)
      colnames(DEG)[1] <- 'gene'
      colnames(DEG)[3] <- 'p_val'
      colnames(DEG)[5] <- 'p_val_adj'
      DEGlist[[paste0('Monocle3_', i)]] <- DEG
    }
  }
  if ('SCDE' %in% test.method) {
    message('##################       Perform SCDE        ##################')
    DEG <- DEG_SCDE(object, groupname = groupname, #batchname  = covariates, 
                    length.out = 400, adjust.method = adjust.method,
                    n.randomizations  =  100,
                    n.cores = n.cores, verbose = 0, ...)
    DEG$gene <- rownames(DEG)
    # colnames(DEG)[4] <- 'avg_log2FC'
    DEGlist[['SCDE']] <- DEG
  }
  if ('miloDE' %in% test.method) {
    message('##################       Perform miloDE        ##################')
    DEG <- DEG_miloDE(object, 
                      groupname = groupname, 
                      covariates = covariates,
                      min.cell = min.cell, 
                      n.cores = n.cores, ...)
    #comparison <- paste(c(ident.2, ident.1), collapse = "_vs_")
    #DEG$de_results$cluster <- comparison
    #colnames(DEG$de_results)[3] <- 'avg_log2FC'
    colnames(DEG$de_results)[4] <- 'p_val'
    colnames(DEG$de_results)[5] <- 'p_val_adj'
    DEGlist[['miloDE']] <- DEG
  }
  if (grepl('DAseq',test.method)) {
    DAseq.method <- gsub('DAseq_','',test.method[grepl('DAseq',test.method)])
    for (i in DAseq.method) {
      message(paste0('##################  Perform DAseq ',i,'  #################'))
      DEG <- DEG_DAseq(object, groupname = groupname, p.adjust.method = adjust.method,
                       ident.1 = ident.1, ident.2 = ident.2, min.cell = min.cell, reduction = 'pca',...)
      DEG$de_results$gene <- rownames(DEG$de_results)
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
                      min.cell = min.cell, 
                      test.meld = test.meld, tmpdir = tmpdir, ...
      )
      DEGlist[[paste0('Meld_',i)]] <- DEG
    }
  }
  
  message('##################  Trimming DEG result...  ##################')
  for (x in names(DEGlist)) {
    if (grepl('DAseq|Meld|milo',x)) {
      degdf <- as.data.frame(DEGlist[[x]]$de_results)
      if (dim( degdf )[1] == 0) next
      meta <- genemeta[unique(degdf$gene),]
      meta$gene <- rownames(meta)
      degdf$avg_log2FC <- NULL
      degdf$pct.1 <- NULL
      degdf$pct.2 <- NULL
      degdf <- left_join(degdf, meta, by = 'gene')
      degdf$cluster <- paste0(ident.1, "_vs_", ident.2)
      DEGlist[[x]]$de_results <- degdf
    }else{
      degdf <- as.data.frame(DEGlist[[x]])
      if (dim( degdf )[1] == 0) next
      meta <- genemeta[unique(degdf$gene),]
      meta$gene <- rownames(meta)
      degdf$avg_log2FC <- NULL
      degdf$pct.1 <- NULL
      degdf$pct.2 <- NULL
      degdf <- left_join(degdf, meta, by = 'gene')
      degdf$cluster <- paste0(ident.1, "_vs_", ident.2)
      DEGlist[[x]] <- degdf
    }
  }
  message('########   Processing ',paste(c(ident.1, test.method), collapse = ', '),' done.  #########')
  return(DEGlist)
}