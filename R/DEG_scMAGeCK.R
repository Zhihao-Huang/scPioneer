#' DEG using scMAGeCK
#' 
#' @param module RRA: test the association with a known marker using scMAGeCK-RRA; LR module for large-scale association test.
#' @param GENE_FRAC A paramter for filtering low expressed genes. By default, only genes that have expressions in at least that fractions of cells (in raw count table) are kept. If raw count table is not found, will look into scaled expression, and only keep genes whose expression is greater than zero in at least that fraction of cells. Default: 0.01 description
#' @examples
#' hmoi <- readRDS('../data/obj_pct01_logfc01.rds')
#' result <- DEG_scMAGeCK(hmoi, covariates = 'Sample')
#' 
#' @export
DEG_scMAGeCK <- function(object, assay_gdo = 'GDO',assay_scripr = 'scripr_thred', 
                         covariates = NULL,
                         BARCODE='/media/london_B/huangzh/project/benchmark_perturb/TEMP/BARCODE_10x.txt',
                         module = c('RRA','LR'),
                         NEGCTRL='non-targeting',target_gene='CD52',
                         label='LR_result' ,PERMUTATION = 100,GENE_FRAC = 0.4){
  module <- match.arg(arg = NULL, choices = module) 
  ## Generate barcode information
  bc_frame <- guidematrix_to_triplet(object@assays[[assay_gdo]]@counts, object)
  meta <- object@assays[[assay_gdo]]@meta.features
  bc_frame[,'gene'] <- meta[bc_frame[,'barcode'], 'grna_target']
  bc_frame[,'sgrna'] <- meta[bc_frame[,'barcode'], 'sgRNA guide sequence']
  ## filter sgRNA
  bn_frame <- guidematrix_to_triplet(object@assays[[assay_scripr]]@counts, object)
  bn_frame <- bn_frame[,1:3]
  colnames(bn_frame)[c(2,3)] <- c('SNP','keep')
  bc_frame$SNP <- gsub('-[0-9]$|-[0-9][0-9]$', '', bc_frame$barcode)
  bc_frame <- left_join(bc_frame, bn_frame)
  sum(is.na(bc_frame$keep))
  bc_frame <- bc_frame[!is.na(bc_frame$keep),]
  bc_frame$SNP <- NULL
  bc_frame$keep <- NULL
  
  write.table(bc_frame,file=BARCODE,row.names = F,sep='\t',quote=F)
  
  feat_c=pre_processRDS(BARCODE = BARCODE, RDS = object)
  feat_c <- NormalizeData(feat_c, normalization.method = "LogNormalize", scale.factor = 10000)
  feat_c <- FindVariableFeatures(feat_c, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(feat_c)
  feat_c <- ScaleData(feat_c, features = all.genes, vars.to.regress = covariates)
  
  #  test the association with a known marker using scMAGeCK-RRA
  if (module == 'RRA') {
    result<-scmageck_rra(BARCODE=bc_frame,RDS=feat_c,GENE=target_gene,
                         LABEL=label,KEEPTMP = FALSE,NEGCTRL=NEGCTRL)
  }
  if (module == 'LR') {
    # LR module for large-scale association test.
    result<-scmageck_lr(BARCODE=BARCODE,RDS=feat_c,LABEL=label,NEGCTRL=NEGCTRL,
                        PERMUTATION = PERMUTATION,GENE_FRAC = GENE_FRAC)
  }
  return(result)
}