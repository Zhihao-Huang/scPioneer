#' DEG using FR-perturb
#' 
#' @param num.perms specifies the number of permutations to perform for permutation testing. Default 1000. 
#' @examples
#' hmoi <- readRDS('../data/obj_pct01_logfc01.rds')
#' DEG_FRperturb(hmoi, assay_crispr = 'crispr_thred', workdir = "/media/london_B/huangzh/project/benchmark_perturb/TEMP/", control.prefix = 'NTC', covariates = c('percent.mt', 'Sample', 'number_positive_sgRNA'))
#' 
#' @export
DEG_FRperturb <- function(object, assay_crispr = 'crispr_thred', workdir = './',
                          run_FR_Perturb_script = '/media/london_B/huangzh/software/FR-Perturb/run_FR_Perturb.py',
                          python2use = "/media/london_B/huangzh/software/miniconda3/envs/FR-Perturb/bin/python",
                          num.perms = NULL,
                          covariates = NULL, control.prefix = NULL) {
  reticulate::use_python(python2use)
  # save Seurat object as anndata h5ad.
  object@assays$RNA@scale.data <- as.matrix(object@assays$RNA@counts[1,,drop = F])
  object@assays$RNA@data <- object@assays$RNA@counts
  SCRNA_h5ad_PATH = paste0(workdir,'sc_RNA.h5ad')
  SCSG_txt_PATH = paste0(workdir,'sc_sg.txt')
  sceasy::convertFormat(object, from="seurat", to="anndata", outFile=SCRNA_h5ad_PATH)
  # save sgRNA binary matrix
  sgmat <- object@assays[[assay_crispr]]$counts
  sgdf <- cbind(data.frame(Cell_Barcode = colnames(sgmat), t(sgmat)))
  write.table(sgdf, quote = F, file = SCSG_txt_PATH, row.names =F)
  
  # run FRperturb on command
  command <- paste0(python2use,' ',run_FR_Perturb_script, ' --input-h5ad ', SCRNA_h5ad_PATH, 
                    ' --input-perturbation-matrix ', SCSG_txt_PATH, ' --compute-pval --out ',
                    workdir)
  if (!is.null(covariates)) command <- paste0(command, ' --covariates ', paste(covariates, collapse = ','))
  if (!is.null(control.prefix)) {
    controls <- colnames(sgdf)[grepl(paste0('^',control.prefix), colnames(sgdf))]
    command <- paste0(command, ' --control-perturbation-name ', paste(controls, collapse = ','))
  }
  if (!is.null(num.perms)) command <- paste0(command, ' --num-perms ', num.perms)
  system(command)
}