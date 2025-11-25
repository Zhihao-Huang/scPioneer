#' DEG using GSFA
#' 
#' @param K Number of factors to use in the model
#' @param prior_type Type of sparse prior used on gene weights, 
#' can be "mixture_normal" or "spike_slab", "mixture_normal" sometimes works better in inducing sparsity.
#' @param init.method Method to initialize the factors, can be one of "svd" (truncated SVD on ‘Y’) or "random"; description
#' @param ... Parameter to GetAssayData.
#' @param G_mat A binary matrix assigning sgRNA for each cell. Each cell contains a unique gRNA. 1 or 0. Column is cell.
#' @examples
#' ## source conda environment GSFA
#' # source activate GSFA 
#' obj <- readRDS('../data/obj_pct01_logfc01.rds')
#' G_mat <- t(as.matrix(obj@@assays$crispr_thred@@counts))
#' obj$Sample <- factor(obj$Sample, levels = unique(obj$Sample))
#' fitdata <- DEG_GSFA(obj, gmat, batchname = 'Sample')
#' 
#' @export
DEG_GSFA <- function(object, G_mat, covariates, 
                     num_top_genes = 6000, 
                     K = 20, prior_type = "mixture_normal", 
                     init.method = "svd",
                     niter = 3000, used_niter = 1000, seed = 14314,...) {
  # 30 mins
  dev_res <- deviance_residual_transform(t(as.matrix(GetAssayData(object, slot = 'counts',...))))
  
  top_gene_index <- select_top_devres_genes(dev_res, num_top_genes = num_top_genes)
  dev_res_filtered <- dev_res[, top_gene_index]
  
  covariate_df <- data.frame(object@meta.data[,covariates],
                             lib_size = object$nCount_RNA,
                             umi_count = object$nFeature_RNA,
                             percent_mt = object$percent.mt)
  dev_res_corrected <- covariate_removal(dev_res_filtered, covariate_df)
  
  scaled.gene_exp <- scale(dev_res_corrected)
  
  sample_names <- colnames(GetAssayData(object, slot = 'counts',...))
  gene_names <- rownames(GetAssayData(object, slot = 'counts',...))
  rownames(scaled.gene_exp) <- sample_names
  colnames(scaled.gene_exp) <- gene_names[top_gene_index]
  
  G_mat <- as.matrix(G_mat)
  # Plot cell number
  num_cells <- colSums(G_mat)
  num_cells_df <- data.frame(locus = names(num_cells),
                             count = num_cells)
  p <- ggplot(data = num_cells_df, aes(x=locus, y=count)) +
    geom_bar(stat = "identity", width = 0.6) +
    labs(x = "Perturbation",
         y = "Number of cells") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size =11),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          panel.grid.minor = element_blank())
  # 30GB, 3 hours
  set.seed(seed)
  fit <- fit_gsfa_multivar(Y = scaled.gene_exp, G = G_mat, 
                           K = K,
                           prior_type = prior_type, 
                           init.method = init.method,
                           niter = niter, used_niter = used_niter,
                           verbose = T, return_samples = T)
  return(fit)
  
}