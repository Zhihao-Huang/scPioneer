#' Extract cells by proportion
#' 
#' @examples 
#' pbmc$Cellname <- colnames(pbmc)
#' subsetID(pbmc@@meta.data[,c('Annotation','Cellname')], num = 100)
#' @export
subsetID <- function (IDmat, num = NULL, fraction = NULL, min.cell = 10, 
                      seed = 123) 
{
  newI <- data.frame(ID = as.vector(IDmat[, 1]), Cellnames = as.vector(IDmat[,2]), 
                     stringsAsFactors = F)
  colnames(newI) <- c("ID", "Cellnames")
  
  if (!is.null(fraction)) {
    set.seed(seed)
    cell10000 <- newI %>% group_by(ID) %>% dplyr::sample_frac(fraction)
    cellnum <-  newI %>% group_by(ID) %>% dplyr::summarise(count = n())
    pos <- cellnum < min.cell
    if (any(pos)) {
      smalltypes <- cellnum[pos,]$ID
      message('The cell number of filtered types ', paste0(paste(smalltypes,collapse = ','), ' is less than min.cell = ', min.cell, '. Keep original cells.'))
      cell10000 <- cell10000[!cell10000$ID %in% smalltypes, ]
      smalldf <- newI[newI$ID %in% smalltypes, ]
      set.seed(seed)
      smallID <- smalldf %>% group_by(ID) %>% dplyr::slice_sample(n = min.cell)
      cell10000 <- rbind(smalldf, smallID)
    }
    if (length(unique(cell10000$ID)) != length(unique(newI$ID))) {
      stop("Cell type lost after subsetting.")
    }
  }else{
    if (is.null(num)) stop('Neither num nor fraction is not default.')
    mingroup <- min(table(newI$ID))
    frac <- round(num/nrow(newI), 5)
    set.seed(seed)
    if (frac > 1) {
      cell10000 <- newI
    }
    else if (mingroup * frac < 1) {
      cell10000 <- lapply(unique(newI$ID), function(x) {
        if (sum(newI$ID %in% x) < 10) {
          newI[newI$ID %in% x, ]
        }
        else if (sum(sample_frac(newI, frac)$ID %in% x) < 
                 min.cell) {
          newI[newI$ID %in% x, ] %>% dplyr::slice_sample(n = min.cell)
        }
        else {
          newI[newI$ID %in% x, ] %>% dplyr::sample_frac(frac)
        }
      })
      cell10000 <- do.call(rbind, cell10000)
    }
    else {
      cell10000 <- newI %>% group_by(ID) %>% dplyr::sample_frac(frac)
    }
  }
  return(cell10000)
}

#' Extract cells by expected number
#' 
#' @examples 
#' subsetID2(newI, expected.cell = 200, seed = 123)
#' @export
subsetID2 <- function (newI, expected.cell = 200, fraction = NULL, min.cell = 10, seed = 123) 
{
  colnames(newI) <- c("ID", "Cellnames")
  if (!is.null(fraction)) {
    cellnum <-  newI %>% group_by(ID) %>% dplyr::summarise(count = n())
    # cell count less than min.cell after fraction subset
    pos <- cellnum$count*fraction < min.cell
    if (any(pos)) {
      smalltypes <- cellnum[pos,]$ID
      message('Warning: the cell number of filtered types ', paste0(paste(smalltypes,collapse = ','), ' is less than min.cell = ', min.cell, '. Restore ', min.cell, ' cells.'))
      cellext <- newI[!newI$ID %in% smalltypes, ]
      set.seed(seed)
      cellext <- cellext %>% group_by(ID) %>% dplyr::sample_frac(fraction)
      smalldf <- newI[newI$ID %in% smalltypes, ]
      set.seed(seed)
      smallID <- smalldf %>% group_by(ID) %>% dplyr::slice_sample(n = min.cell)
      cellext <- rbind(cellext, smallID)
    }else{
      set.seed(seed)
      cellext <- newI %>% group_by(ID) %>% dplyr::sample_frac(fraction)
    }
    if (length(unique(cellext$ID)) != length(unique(newI$ID))) {
      stop("Cell type lost after subsetting.")
    }
  }else{
    mingroup <- min(table(newI$ID))
    maxgroup <- max(table(newI$ID))
    if (mingroup > expected.cell) {
      set.seed(seed)
      cellext <- newI %>% group_by(ID) %>% dplyr::slice_sample(n = expected.cell)
    }
    else if (maxgroup < expected.cell) {
      message(paste0("All the cell numbers of cluster are less than ", 
                     expected.cell))
      cellext <- newI
    }
    else {
      smallgroupnames <- names(table(newI$ID))[table(newI$ID) < 
                                                 expected.cell]
      smallgroup <- newI[newI$ID %in% smallgroupnames, ]
      largegroup <- newI[!newI$ID %in% smallgroupnames, ]
      set.seed(seed)
      largegroupext <- largegroup %>% group_by(ID) %>% dplyr::slice_sample(n = expected.cell)
      cellext <- rbind(smallgroup, as.matrix(largegroupext))
    }
  }
  return(cellext)
}


