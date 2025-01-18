#' Get genelist with symbol and ensembl ID.
#'
#' @param version Version of ensembl.
#' @param using.local.file Whether to use annotation file from local directory instead of connecting to internet.
#' @param file.dir Directory for annotation file.
#' @param remove.duplicated Whether to remove duplicated gene symbol.
#' @param only.one.to.one Whether to keep rows in which both gene symbol and ensembl ID are unique.
#'
#' @examples
#' df <- Get_symbol_ensembl_ID(version = 98)
#' df <- Get_orthologs_mouse_human(version = 98, using.local.data = T) # if cannot access biomaRt.
#'
#' @export
Get_symbol_ensembl_ID <- function (version = 98, species = "Human", 
                                   using.local.file = F, file.dir = NULL, 
                                   using.local.data = F,
                                   remove.duplicated = T, only.one.to.one = F) 
{
  if (using.local.file) {
    gene_anno <- read.table(file = paste0(file.dir, "symbol_ensembl_ID_version_", 
                                          version, ".txt"), sep = "\t")
  }else if (using.local.data) {
    gene_anno <- symbol_ensembl_list[[species]][[paste0('v',version)]]
  }else {
    if (species == "Human") {
      ensembl <- tryCatch(biomaRt::useEnsembl(biomart = "ensembl", 
                                              dataset = "hsapiens_gene_ensembl", version = version), 
                          error = function(cond) {
                            message("Could not connect to the ensembl URL, recommend to use local file by setting using.local.file = T.")
                            message("Here's the original error message:")
                            stop(cond)
                          })
    }
    else if (species == "Mouse") {
      ensembl <- tryCatch(biomaRt::useEnsembl(biomart = "ensembl", 
                                              dataset = "mmusculus_gene_ensembl", version = version), 
                          error = function(cond) {
                            message("Could not connect to the ensembl URL, recommend to use local file by setting using.local.file = T.")
                            message("Here's the original error message:")
                            stop(cond)
                          })
    }
    else {
      stop("Only human or mouse is allowed.")
    }
    gene_anno <- biomaRt::getBM(attributes = c("ensembl_gene_id", 
                                               "hgnc_symbol", "chromosome_name", "start_position", 
                                               "end_position"), filters = "chromosome_name", values = c(as.character(1:22), 
                                                                                                        "X", "Y"), mart = ensembl)
    colnames(gene_anno) <- c("ensembl_gene_id", "symbol", 
                             "chromosome_name", "start_position", "end_position")
  }
  if (remove.duplicated) {
    pos_empty <- gene_anno$symbol != ""
    gene_anno_filter <- gene_anno[pos_empty, ]
    pos_redup <- !duplicated(gene_anno_filter$ensembl_gene_id) & 
      !duplicated(gene_anno_filter$symbol)
    gene_anno_filter <- gene_anno_filter[pos_redup, ]
  }
  else if (only.one.to.one) {
    duplicated_symbol <- unique(gene_anno$symbol[duplicated(gene_anno$symbol)])
    duplicated_ensembl <- unique(gene_anno$ensembl_gene_id[duplicated(gene_anno$ensembl_gene_id)])
    gene_anno_filter <- gene_anno[!gene_anno$symbol %in% 
                                    duplicated_symbol & !gene_anno$ensembl_gene_id %in% 
                                    duplicated_ensembl, ]
  }
  else {
    gene_anno_filter <- gene_anno
  }
  return(gene_anno_filter)
}
