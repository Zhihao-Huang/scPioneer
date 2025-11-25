#' Check whether there is special symbol in names of cell barcode or annotation to avoid error.
#'
#' @param celltype A vector, can be barcode or annotation.
#'
#' @export
check_replace_name <- function(celltype, pattern = NULL, replacement = NULL, fixed = F) {
  # user pattern
  if (!is.null(pattern) & !is.null(replacement)) {
    if (any(grepl(pattern = pattern, x = celltype, fixed = fixed))) {
      warning(paste0("Feature names cannot have characters (",pattern,"), replacing with ampersand (",replacement,")"),
              call. = FALSE, immediate. = TRUE)
      celltype <- gsub(pattern = pattern, replacement = replacement,
                       x = celltype, fixed = fixed)
      
    }
  }
  return(celltype)
}

#' Filter genes by min.pct and logfc.threshold
#' 
#' @references Seurat::FoldChange
#' @export
filter_gene <- function(object, group.by = NULL, ident.1, ident.2 = NULL, 
                        min.pct = 0.01, logfc.threshold = 0,
                        assay = 'RNA', slot = 'counts', return.fcresult = F) {
  # filtering genes
  fcresult <- Seurat::FoldChange(object, ident.1, 
                                 group.by = group.by, assay = assay, slot = slot)
  # feature selection (based on percentages)
  alpha.min <- pmax(fcresult$pct.1, fcresult$pct.2)
  names(x = alpha.min) <- rownames(x = fcresult)
  features <- names(x = which(x = alpha.min >= min.pct))
  fcresult <- fcresult[features,]
  fcresult <- fcresult[abs(fcresult$avg_log2FC) >= logfc.threshold, ]
  message(paste0('Number of genes remained: ', nrow(fcresult)))
  if (return.fcresult) {
    return(fcresult)
  }else{
    object <- object[rownames(fcresult),]
    return(object)
  }
}