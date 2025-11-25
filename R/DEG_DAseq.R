#' Pairwise DEG using DAseq
#' 
#' @param batchname Colname including covariate in meta.data.
#' @param groupname Colname including two groups for DEG analysis, e.g. SNP and WT labels.
#' @param k.vector k values to create the score vector.
#' @param k.folds integer, number of data splits used in the neural network, default 10.
#' @param n.runs integer, number of times to run the neural network to get the predictions, default 5.
#' @param n.rand integer, number of random permutations to run, default 2.
#' @param pred.thres.range length-2 vector, top and bottom threshold on DA measure, default NULL, select significant DA cells based on permutation
#' @param deg.func STG or Seurat. STG: Use STG (stochastic gates) to select genes that separate each DA region from the rest of the cells.
#' @param lambda numeric, regularization parameter that weights the number of selected genes, a larger lambda leads to fewer genes, default 1.5
#' @param ... Parameter to GetAssayData.
#' 
#' @examples
#' packageVersion('Seurat')
#' # ‘4.4.0’
#' packageVersion('Matrix')
#' # '1.6.2'
#' pbmc <- readRDS('../data/pbmc_v4.rds')
#' pbmc <- subset(pbmc, Annotation %in% c('B','NK'))
#' pbmc$Annotation <- factor(pbmc$Annotation, levels = unique(pbmc$Annotation))
#' # includes random effects
#' pbmc$batch <- c(rep('sample1',100), rep('sample2',100),rep('sample3',299))
#' DEG <- DEG_DAseq(pbmc, groupname = 'Annotation', ident.1 = 'B', ident.2 = 'NK')
#' 
#' @export
DEG_DAseq <- function(object, groupname, ident.1, ident.2,
                      lambda = 1.5, n.cores = 4, return.model = T,
                      da.regions.to.run = NULL,
                      p.adjust.method = 'fdr',
                      min.pct = 0.2, logfc.threshold = 0,
                      min.cell = 50,
                      python2use = "/media/london_B/huangzh/software/miniconda3/envs/FR-Perturb/bin/python",
                      GPU = 4, reduction = 'pca',...) {
  # filter genes
  fc_results <- filter_gene(object = object, group.by = groupname,
                            ident.1 = ident.1, ident.2 = ident.2, min.pct = min.pct, 
                            logfc.threshold = logfc.threshold, slot = 'data', 
                            return.fcresult = T)
  fc_results$gene <- rownames(fc_results)
  object <- object[rownames(fc_results),]
  
  X.label.objectdata <- as.character(object@meta.data[,groupname])
  object$X.label.objectdata <- X.label.objectdata
  labels_res <- ident.1
  labels_nonres <-ident.2
  
  object <- subset(object, X.label.objectdata %in% c(labels_res, labels_nonres))
  X.data.objectdata <- GetAssayData(object, slot = 'data',...)
  X.objectdata <- object@reductions$pca@cell.embeddings
  X.2d.objectdata <- object@reductions[[reduction]]@cell.embeddings[,1:2] # only for purpose of visualization
  X.label.objectdata <- object$X.label.objectdata
  
  da_cells <- getDAcells(
    X = X.objectdata,
    cell.labels = X.label.objectdata,
    labels.1 = labels_res,
    labels.2 = labels_nonres,
    plot.embedding = X.2d.objectdata
  )
  if (length(da_cells$da.up) < min.cell | length(da_cells$da.down) < min.cell ) {
    message(paste0("DA cells in up or down group has cells fewer than ",min.cell," cells. Return empty."))
    return(list(object = list(da_cells = da_cells),
                de_results = NULL))
  }
  da_regions <- getDAregion(
    X = X.objectdata,
    da.cells = da_cells,
    cell.labels = X.label.objectdata,
    labels.1 = labels_res,
    labels.2 = labels_nonres,
    plot.embedding = X.2d.objectdata, 
    do.plot = F
  )
  tab <- table(da_regions$da.region.label)
  exist_DA <- F
  if (length(tab) > 2) exist_DA <- tab[2] > 2 | tab[3] > 2
  if (!exist_DA) {
    message(" Up or down DA regions has cells fewer than 3 cells. Retrun empty.")
    return(list(object = list(da_cells = da_cells,  da_regions = da_regions),
                de_results = NULL))
  }

  # t.test and wilcox test
  assignInNamespace("is_conda_python", function(x){ return(FALSE) }, ns="reticulate")
  STG_markers <- STGmarkerFinder(
    X = X.data.objectdata,
    da.regions.to.run = da.regions.to.run,
    da.regions = da_regions,
    lambda = lambda, n.runs = n.cores, return.model = return.model,
    python.use = python2use, GPU = GPU
  )
  deglist <- STG_markers$da.markers
  deglist <- lapply(deglist, function(x) {
    colnames(x)[3]  <- 'p_val'
    x
  })
  
  DEG <- do.call(rbind, deglist)
  DEG <- as.data.frame(DEG)
  message(paste0(length(unique(DEG$DA.region)), ' DA regions left.'))
  DEG$p_val_adj <- p.adjust(DEG$p_val, method = p.adjust.method)
  DEG <- left_join(DEG[, c('p_val','p_val_adj', 'gene')], fc_results, by = 'gene')
  DEG <- DEG[,c('p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj','gene')]
  DEG <- DEG %>% mutate(cluster = ifelse(avg_log2FC < 0, 2, 1)) %>% arrange(cluster,p_val_adj) %>% mutate(cluster = ifelse(cluster == 1, ident.1, ident.2))
  results <- list(object = list(da_cells = da_cells, da_regions = da_regions),
                  de_results = DEG)
  return(results)
}
