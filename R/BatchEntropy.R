#' Calculate the cell entropy according to the cell group. (From Liziyi, onco-fetal_utils.R)
#'
#' @param input.data A matrix with each cell in a row. The column could be
#' tSNE coordinate, expression matrix or even a dist format.
#' @param group.id A vector with the same length as the row of input.data, or
#' a list of vectors.
#' @param k.used The k used to build a kNN-graph.
#' @param dimension.used The method to reduce the dimension, tSNE by default,
#' could also be PCA or raw.
#'
#' @author Liziyi
#' @return A vector with each cell's entropy.
#' @export

BatchEntropy <- function(input.data, group.id, 
                         k.used = 30, 
                         dimension.used = "raw",
                         n.cores = 1) {
  library(dbscan)
  library(entropy)
  library(parallel)
  if (dimension.used == "raw") {
    knn_graph <- kNN(x = input.data, k = k.used, sort = FALSE, search = "dist")
  }
  if (dimension.used == "tSNE") {
    cat(paste("Reduce the dimension using", dimension.used, "at", Sys.time(), "\n"))
    tSNE.coor <- Rtsne::Rtsne(input.data)
    cat(paste("Buld a k-NN graph at", Sys.time(), "\n"))
    knn_graph <- kNN(x = tSNE.coor$Y, k = k.used, sort = FALSE, search = "dist")
  }
  if (dimension.used == "PCA") {
    cat(paste("Reduce the dimension using", dimension.used, "at", Sys.time(), "\n"))
    PCA.coor <- prcomp(input.data, center = FALSE)
    PCA.cumsd <- cumsum(PCA.coor$sdev) / sum(PCA.coor$sdev)
    nPCs.used <- which(PCA.cumsd > 0.9)[1]
    cat(paste("Buld a k-NN graph at", Sys.time(), "\n"))
    knn_graph <-
      kNN(x = PCA.coor$x[, 1:nPCs.used], k = k.used, sort = FALSE, search = "dist")
  }
  if (!is.list(group.id)) {
    group.id <- list(default = group.id)
  }
  cell_entropy <- mclapply(names(group.id), function(i) {
    knn_group <- matrix(group.id[[i]][knn_graph$id],
                        nrow = nrow(input.data),
                        byrow = FALSE)
    row.names(knn_group) <- row.names(input.data)
    colnames(knn_group) <- 1:k.used
    cat(paste("Calculate the cell entropy of", i, "at", Sys.time(), "\n"))
    apply(knn_group, 1, function(x) {entropy(table(x))})
  }, mc.cores = n.cores)
  names(cell_entropy) <- names(group.id)
  return(cell_entropy)
}
