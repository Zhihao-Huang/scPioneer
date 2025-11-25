#' Pairwise DEG using Monocle3
#' 
#' @param batchname Colname including covariate in meta.data.
#' @param groupname Colname including two groups for DEG analysis, e.g. SNP and WT labels.
#' @param expression_family DEG models passed to fit_models().
#' @param ... Parameter to GetAssayData.
#' @examples
#' pbmc <- readRDS('../data/pbmc_v4.rds')
#' pbmc <- subset(pbmc, Annotation %in% c('B','NK'))
#' pbmc$Annotation <- factor(pbmc$Annotation, levels = unique(pbmc$Annotation))
#' # includes random effects
#' pbmc$batch <- c(rep('sample1',100), rep('sample2',100),rep('sample3',299))
#' DEG <- DEG_Monocle3(pbmc, groupname = 'Annotation',  batchname  = 'batch')
#' 
#' @export
DEG_Monocle3 <- function(object,
                         groupname, ident.1 = ident.1, ident.2 = ident.2,
                         batchname = NULL,
                         min.pct = 0.2, logfc.threshold = 0,
                         expression_family = c("quasipoisson","negbinomial",'poisson','binomial',
                                               "gaussian", "zipoisson", "zinegbinomial", "mixed-negbinomial"),
                         adjust.method = c('bonferroni','holm', 'hochberg', 'hommel', 'BH', 'BY', 'fdr', 'none'),
                         n.cores = 8,
                         ...) {
  adjust.method <- match.arg(arg = NULL, choices = adjust.method)
  expression_family <- match.arg(arg = NULL, choices = expression_family)
  cellID <- colnames(object)[object@meta.data[,groupname] %in% c(ident.1, ident.2)]
  object <- subset(object, cells = cellID)
  # filter genes
  fc_results <- filter_gene(object = object, group.by = groupname,
                            ident.1 = ident.1, ident.2 = ident.2, min.pct = min.pct, 
                            logfc.threshold = logfc.threshold, slot = 'data', 
                            return.fcresult = T)
  fc_results$gene <- rownames(fc_results)
  object <- object[rownames(fc_results),]
  expression_matrix <- GetAssayData(object, slot = 'counts',...)
  cell_metadata <- object@meta.data
  gene_annotation <- data.frame(id = rownames(object), gene_short_name = rownames(object), 
                                num_cells_expressed = apply(expression_matrix, 1, function(x) sum(x >0)))
  cds <- new_cell_data_set(expression_matrix, cell_metadata = cell_metadata, 
                           gene_metadata = gene_annotation)
  if (is.null(batchname)) {
    model_formula_str = paste0('~',groupname)
  }else{
    model_formula_str = paste0('~',groupname,' + ', paste(batchname, collapse = ' + '))
  }
  gene_fits <- fit_models(cds, expression_family = expression_family,
                          cores = n.cores, model_formula_str = model_formula_str)
  fit_coefs <- coefficient_table(gene_fits)
  fit_coefs = fit_coefs %>% dplyr::group_by(model_component, term) %>% 
    dplyr::mutate(adj.P.Val = stats::p.adjust(p_value, method = adjust.method)) %>% 
    dplyr::ungroup()
  print(head(fit_coefs[,8:16]))
  deg <- fit_coefs %>% filter(term != "(Intercept)") %>%
    select(gene_short_name, term, p_value, q_value, adj.P.Val, estimate)
  deg <- deg[grepl(groupname,deg$term),]
  colnames(deg)[1] <- 'gene'
  colnames(deg)[3] <- 'p_val'
  colnames(deg)[5] <- 'p_val_adj'
  print(head(deg))
  DEG <- left_join(deg[, c('p_val','p_val_adj', 'gene')], fc_results, by = 'gene')
  DEG <- DEG[,c('p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj','gene')]
  DEG <- DEG %>% mutate(cluster = ifelse(avg_log2FC < 0, 2, 1)) %>% arrange(cluster,p_val_adj) %>% mutate(cluster = ifelse(cluster == 1, ident.1, ident.2))
  return(DEG)
}