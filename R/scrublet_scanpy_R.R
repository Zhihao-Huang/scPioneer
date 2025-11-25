#' @title Scrublet
#' @description See article doi: https://doi.org/10.1016/j.cels.2018.11.005
#' @param seurat_obj the CellDataSet upon which to perform Scrublet
#' @param outdir Derectory to save the result and temp result of scrublet
#' @param python_home The python home directory where Scrublet is installed
#' @param return_results_only bool (optional, default False)
#' @param source_py_script scrublet python script. default is script in inst/
#' @param batch_key colname including sample names for batch correlation.  
#' @param expected_doublet_rate, float (optional, default=0.05), See scrublet reference - expected_doublet_rate: the 
#' fraction of transcriptomes that are doublets, typically 0.05-0.1. Results are not particularly sensitive to this parameter. For this example, the expected doublet rate comes from the Chromium User Guide: https://support.10xgenomics.com/permalink/3vzDu3zQjY0o2AqkkkI4CC

#' @return A data.frame including Cellname,doublet_score,predicted_doublet,Sample
#' @importFrom reticulate use_python 
#' @importFrom reticulate source_python
#' @examples 
#' pbmc$Sample <- 'single'
#' res <- scrublet_scanpy_R(pbmc, outdir = './test/')
#' 
#' @export
scrublet_scanpy_R <- function (seurat_obj, outdir,
                               python_home = system("which python", intern = TRUE), 
                               source_py_script = NULL,
                               expected_doublet_rate=0.05, batch_key = 'None',
                               return_results_only = T) 
{
  print(reticulate::py_config())
  if(!reticulate::py_module_available("scanpy")) {
    message("Python package scanpy is not be installed; - try running 'py_config()'")
    message("install scanpy by py_install('scanpy')...")
    reticulate::py_install('scanpy')
  }
  if(!reticulate::py_module_available("scrublet")) {
    message("Python package scrublet is not be installed; - try running 'py_config()'")
    message("install scrublet by py_install('scrublet')...")
    reticulate::py_install('scrublet')
  }
  if (is.null(source_py_script)) {
    source_py_script <- system.file("scrublet_scanpy.py", package = "scPioneer")
  }
 
  reticulate::source_python(source_py_script)
  
  adata_path = paste0(outdir,'adata_scrublet.h5ad')
  df_path = paste0(outdir,'df_scrublet.csv')
  
  rds_to_h5ad(seurat_obj, h5ad_path = adata_path)
  scrublet_scanpy_args<-c(list(adata_path, out_path = df_path, batch_key = batch_key,
                               expected_doublet_rate = expected_doublet_rate))
  do.call(scrublet_scanpy, scrublet_scanpy_args)
  scrublet_res <- read.csv(df_path)
  scrublet_res <- scrublet_res[, c('Cellname','doublet_score','predicted_doublet','Sample')]
  names(scrublet_res)<-c("doublet_scores", "predicted_doublets")
  if (return_results_only) {
    return(scrublet_res)
  }
  else {
    seurat_obj$doublet_scores<-scrublet_res$doublet_scores
    seurat_obj$predicted_doublets<-scrublet_res$predicted_doublets
    return(seurat_obj)
  }
}

rds_to_h5ad <- function(object, h5ad_path, assay = 'RNA') {
  obj_raw <- object
  obj_raw[[assay]] <- as(obj_raw[[assay]], "Assay")
  obj_raw@assays[[assay]]@scale.data <- as.matrix(1)
  obj_raw@assays[[assay]]@data <- GetAssayData(obj_raw, assay = assay, layer = 'counts')
  sceasy::convertFormat(obj_raw, from="seurat", to="anndata", outFile= h5ad_path)
}