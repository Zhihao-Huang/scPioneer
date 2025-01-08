#' Get orthologs between mouse and human.
#'
#' @param version Version of ensembl.
#' @param using.local.file Whether to use annotation file from local directory instead of connecting to internet.
#' @param file.dir Directory for annotation file.
#' @param remove.duplicated Whether to remove duplicated gene symbol.
#' @param only.one.to.one Whether to keep rows in which both gene symbol and ensembl ID are unique.
#'
#' @examples
#' df <- Get_orthologs_mouse_human(using.local.file=T)
#' filedir <- '/BGFS1/home/huangzh/workspace/database/orthologs/mouse_human/'
#' df <- Get_orthologs_mouse_human(version = 104, using.local.file = T,file.dir = filedir, only.one.to.one = T)
#'
#' @export
Get_orthologs_mouse_human <- function(version = 104,
                                      using.local.file = F,
                                      file.dir = '/BGFS1/home/huangzh/workspace/database/orthologs/mouse_human/',
                                      using.local.data = F,
                                      remove.duplicated = T,
                                      only.one.to.one = F) {
  if (using.local.file) {
    geneset <- read.table(file = paste0(file.dir,'orthologs_mouse_human_version_',version,'.txt'),sep = '\t')
  }else if (using.local.data) {
    geneset <- orthologslist[[paste0('v',version)]]
  }else{
    human <- tryCatch(biomaRt::useEnsembl(biomart = 'ensembl',
                                          dataset = 'hsapiens_gene_ensembl',
                                          mirror = 'asia',
                                          version = version),
                      error=function(cond) {
                        message("Could not connect to the ensembl URL, recommend to use local file by setting using.local.file = T.")
                        message("Here's the original error message:")
                        stop(cond)
                        # Choose a return value in case of error
                        #return(NA)
                      })
    mouse <- biomaRt::useEnsembl(biomart = 'ensembl',
                                 dataset = 'mmusculus_gene_ensembl',
                                 mirror = 'asia',
                                 version = version)
    #gene_positions <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position"),
    #                        filters = "chromosome_name", values = 'Y', mart = human)
    geneset <- biomaRt::getLDS(attributes = c("hgnc_symbol","ensembl_gene_id"),
                               mart = human,
                               attributesL = c("mgi_symbol","ensembl_gene_id"),
                               martL = mouse,
                               filters = "chromosome_name",
                               values = c(as.character(1:22),'X','Y'),
                               uniqueRows= F)
    colnames(geneset) <- c('HGNC_symbol','HGNC_ensembl_gene_id','MGI_symbol','MGI_ensembl_gene_id')
  }
  if (remove.duplicated) {
    pos_empty <- geneset$HGNC_symbol != '' & geneset$MGI_symbol != ''
    geneset_filter <- geneset[pos_empty,]
    pos_redup <-  !duplicated(geneset_filter$MGI_symbol) &
      !duplicated(geneset_filter$HGNC_symbol)
    geneset_filter <- geneset_filter[pos_redup,]
  } else if (only.one.to.one) {
    duplicated_MGI_gene <- unique(geneset$MGI_symbol[duplicated(geneset$MGI_symbol)])
    duplicated_HGNC_gene <- unique(geneset$HGNC_symbol[duplicated(geneset$HGNC_symbol)])
    geneset_filter <- geneset[!geneset$MGI_symbol %in% duplicated_MGI_gene &
                                !geneset$HGNC_symbol %in% duplicated_HGNC_gene,]
    geneset_filter <- geneset_filter[!duplicated(geneset_filter$HGNC_ensembl_gene_id),]
  }else{
    geneset_filter <- geneset
  }

  return(geneset_filter)
}

#' Convert symbol gene names between mouse and human.
#' @param geneset Vector of genes.
#' @param species Only human or mouse.
#' @param version Version of ensembl.
#' @param using.local.file Whether to use annotation file from local directory instead of connecting to internet.
#' @param file.dir Directory for annotation file.
#' @param remove.duplicated Whether to remove duplicated gene symbol.
#' @param only.one.to.one Whether to keep rows in which both gene symbol and ensembl ID are unique.
#'
#' @examples
#' df <- Get_orthologs_mouse_human(using.local.file=T)
#' filedir <- '/BGFS1/home/huangzh/workspace/database/orthologs/mouse_human/'
#' df <- Get_orthologs_mouse_human(version = 104, using.local.file = T,file.dir = filedir, only.one.to.one = T)
#'
#' @export
convert_symbol_mouse_human <- function(geneset,
                                       species,
                                       version = 101,
                                       using.local.file = F,
                                       file.dir = '/BGFS1/home/huangzh/workspace/database/orthologs/mouse_human/',
                                       using.local.data = F,
                                       remove.duplicated = T,
                                       only.one.to.one = F
) {
  orthologs <- Get_orthologs_mouse_human(version = version,
                                         file.dir = file.dir,
                                         remove.duplicated = remove.duplicated,
                                         only.one.to.one = only.one.to.one,
                                         using.local.file = using.local.file,
                                         using.local.data = using.local.data)
  if (species == 'mouse') {
    df <- data.frame(MGI_symbol = geneset, stringsAsFactors = F)
    df <- left_join(df, orthologs)
  }else if (species == 'human') {
    df <- data.frame(HGNC_symbol = geneset, stringsAsFactors = F)
    df <- left_join(df, orthologs)
  }
  return(df)
}
