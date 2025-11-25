#' Pairwise DEG using miloDE. Employ edgeR testing (using glmQLFit) within each neighbourhood.
#' Refer to https://github.com/MarioniLab/miloDE?tab=readme-ov-file
#' 
#' @param k defines how many neighbours to use for the neighbourhood assignment. Default k = 25.
#' @param order In c(1,2), defines which order of neighbours to use. Default order = 2.
#' @param ... Parameter to GetAssayData.
#' @examples
#' obj <- readRDS('./data/obj_pct01_logfc01.rds')
#' DEG <- DEG_miloDE(object = obj, samplename = 'Sample',groupname = 'crispr', min.cell = 20,covariates = c('crispr','Sample','percent.mt'))
#' 
#' @export
DEG_miloDE <- function(object, groupname, ident.1, ident.2, covariates = NULL, samplename = NULL, 
                       min.pct = 0.2, logfc.threshold = 0,
                       reducedDim_name = "PCA", k = 25, order = 2, 
                       min_n_cells_per_sample = 3, min_count = 3, min.cell = 50, 
                       filtering = TRUE, n.cores = 1) {
  if (is.null(samplename)) samplename <- 'orig.ident'
  # filter genes
  fc_results <- filter_gene(object = object, group.by = groupname,
                            ident.1 = ident.1, ident.2 = ident.2, min.pct = min.pct, 
                            logfc.threshold = logfc.threshold, slot = 'data', 
                            return.fcresult = T)
  fc_results$gene <- rownames(fc_results)
  object <- object[rownames(fc_results),]
  object$groupname <- object@meta.data[,groupname]
  object$groupname <- factor(object$groupname, levels = c(ident.1, ident.2))
  sce <- Seurat::as.SingleCellExperiment(object)
  # Neighbourhood assignment
  scemilo = assign_neighbourhoods(sce, k = k, order = order, 
                                  filtering = filtering, reducedDim_name = reducedDim_name)
  nhoods_mat <- nhoods(scemilo)
  nhoods_group <- as.matrix(nhoods_mat) * as.numeric(object$groupname)
  nhoods_check <- apply(nhoods_group, 2, function(x) any(table(as.character(x)) < min.cell))
  nhoods_check <- as.vector(nhoods_check)
  if (any(nhoods_check)) {
    warning(paste0('Neighbourhoods ',paste(which(nhoods_check), collapse = ','), ' have groups with cells fewer than ', min.cell, '. Descarded.'))
  }
  # DE testing
  DEG = de_test_neighbourhoods(scemilo, sample_id = samplename, subset_nhoods = which(!nhoods_check),
                               min_n_cells_per_sample = min_n_cells_per_sample,
                               min_count = min_count,
                               design = ~groupname, 
                               covariates = c('groupname', covariates),
                               BPPARAM = MulticoreParam(n.cores),
                               verbose = TRUE)
  DEG <- as.data.frame(DEG)
  colnames(DEG)[4] <- 'p_val'
  colnames(DEG)[5] <- 'p_val_adj'
  DEG <- left_join(DEG[, c('p_val','p_val_adj', 'gene')], fc_results, by = 'gene')
  DEG <- DEG[,c('p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj','gene')]
  DEG <- DEG %>% mutate(cluster = ifelse(avg_log2FC < 0, 2, 1)) %>% arrange(cluster,p_val_adj) %>% mutate(cluster = ifelse(cluster == 1, ident.1, ident.2))
  Milo_results <- list(obj = scemilo, de_results = DEG)
  return(Milo_results)
}