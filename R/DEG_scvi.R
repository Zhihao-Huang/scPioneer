#' Pairwise DEG using scvi
#' 
#' @param on_disk Whether to save temporary results on disk.
#' @param test.Seurat Employ DE methods built in Seurat.
#' @param test.meld Employ DE methods built in Meld.
#' @param return.classes whether to return clusters of cells.
#' @param ... Parameter to GetAssayData.
#' 
#' @references https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/scVI_DE_worm.html
#' @examples
#' pbmc <- readRDS('/media/london_B/huangzh/project/benchmark_perturb/data/pbmc_v4.rds')
#' pbmc <- subset(pbmc, Annotation %in% c('B','NK'))
#' pbmc$Annotation <- factor(pbmc$Annotation, levels = unique(pbmc$Annotation))
#' # includes random effects
#' pbmc$batch <- c(rep('sample1',100), rep('sample2',100),rep('sample3',299))
#' DEG <- DEG_scvi(pbmc, 'Annotation', 'B','NK')
#' DEG <- DEG_scvi(pbmc, 'Annotation', 'B','NK', covariates = c('batch'))
#' 
#' @export
DEG_scvi <- function(object, groupname, ident.1, ident.2,
                     covariates = NULL, n.cores = 8,
                     min.pct = 0.2, logfc.threshold = 0,
                     python2use = '/lustre/home/kwxiong/Huangzhihao/software/miniconda3/envs/scvi-env/bin/python',
                     scvi_script_path = '/lustre/home/kwxiong/Huangzhihao/project/benchmark_perturb//shell/Perturb_benchmark-main/script/scvi.py'
) {
  object <- subset(object, cells = colnames(object)[object@meta.data[,groupname] %in% c(ident.1, ident.2)])
  fc_results <- filter_gene(object = object, group.by = groupname,
                            ident.1 = ident.1, ident.2 = ident.2, min.pct = min.pct, 
                            logfc.threshold = logfc.threshold, slot = 'data', 
                            return.fcresult = T)
  # filter genes
  object <- object[rownames(fc_results),]
  assignInNamespace("is_conda_python", function(x){ return(FALSE) }, ns="reticulate")
  if (!is.null(python2use)) reticulate::use_python(python2use)
  print(reticulate::py_config())
  installed <- reticulate::py_module_available("scvi")
  if(!installed) stop("python module scvi does not seem to be installed; - try running 'py_config()'")
  reticulate::source_python(scvi_script_path)
  #### Input expression matrix
  mat <- GetAssayData(object, slot = 'counts')
  print(class(mat))
  print(dim(mat))
  #X <- bigdataTOmat(mat, min_num = 20000, transpose = T, as.type = "TsparseMatrix")
  X <- as(Matrix::t(mat), "TsparseMatrix")
  i <- as.integer(X@i)
  j <- as.integer(X@j)
  val <- X@x
  dim <- as.integer(X@Dim)
  #### Input meta data
  meta <- object@meta.data[,c(groupname, covariates), drop = F]
  group <- as.vector(meta[,c(groupname)])
  group <- gsub(ident.1,'stimulated', group)
  group <- gsub(ident.2,'unstimulated', group)
  meta$condition <- group
  meta$cell <- rownames(meta)
  if (is.null(covariates)) {
    covariates <- 'None'
  }
  DEG_scvi_py_args<-c(list(i=i, j=j, val=val, dim=dim, cell_idx=colnames(object),
                           feature_idx = rownames(object),
                           meta=meta, batch_key=covariates, n_cores = n.cores))
  
  de_results <- do.call(DEG_scvi_py, DEG_scvi_py_args)
  return(de_results)
}
