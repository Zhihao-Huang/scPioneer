#' Shuffle cells for iSEE
#'
#' @param SeuratS4 A Seurat object.
#' @param num Shuffle cell number.
#' @param select.meta metadata reserved.
#' @return
#' A SingleCellExperiment object.
#' @examples
#' data('pbmc')
#' pbmc@@meta.data$Annotaion <- Idents(pbmc)
#' m_s4 <- sce_10000(SeuratS4 = pbmc, num = 1000, select.meta = c("Annotaion","nFeature_RNA",
#' "nCount_RNA"))
#' @export

sce_10000 <- function(SeuratS4 = SeuratS4, num = 1000, select.meta = c("nFeature_RNA","nCount_RNA"),
                      Assay_used = "RNA",reduction = 'tsne') {
  ## 1. extract < 1000 cells newI is a factor of celltypes
  newI <- Idents(SeuratS4)
  cell_ident <- as.data.frame(newI)
  cell_ident$ID <- rownames(cell_ident)
  mingroup <- min(table(newI))
  frac <- round(num/length(newI), 5)
  if (frac > 1) {
    cell10000 <- cell_ident
  } else if (mingroup * frac < 1) {
    # small group may have cell < 1 after subset, so increase the frac.
    frac <- 1/(mingroup - 1)
    cell10000 <- cell_ident %>% group_by(newI) %>% sample_frac(frac)
  } else {
    cell10000 <- cell_ident %>% group_by(newI) %>% sample_frac(frac)
  }
  
  Idents(SeuratS4) <- factor(Idents(SeuratS4), levels = sort(levels(Idents(SeuratS4))))
  DefaultAssay(SeuratS4) <- Assay_used
  ## 2. subset data and convert it to sce, then save as rds
  m_s4 <- subset(SeuratS4, cells = cell10000$ID)
  m_s4@meta.data <- m_s4@meta.data[, select.meta]
  m_s4@reductions$pca <- NULL
  if (reduction == 'tsne') {
    m_s4@reductions$umap <- NULL
  }else if (reduction == 'umap') {
    m_s4@reductions$tsne <- NULL
  }
  m_sce <- as.SingleCellExperiment(m_s4)
  m_sce@colData$ident <- NULL
  SingleCellExperiment::counts(m_sce) <- SeuratS4@assays[[Assay_used]]@data[, colnames(m_sce)]
  return(m_sce)
}
#' Visualize single cell RNA-Seq data in HTML format.
#'
#' @param data A SingleCellExperiment object.
#' @examples
#' data('pbmc')
#' pbmc@@reductions$pca <- NULL
#' m_sce <- as.SingleCellExperiment(pbmc)
#' SingleCellExperiment::counts(m_sce) <- pbmc@@assays$RNA@@data
#' m_sce@@colData$Annotation <- Idents(pbmc)
#' #iSEEplot(m_sce)
#' @export
iSEEplot <- function(data){
  initialPanels <- DataFrame(
    Name=c("Reduced dimension plot 1", "Reduced dimension plot 2", "Feature assay plot 1"),
    Width=c(4L, 4L, 4L),
    Height=c(600L, 600L, 600L)
  )
  
  rdArgs <- iSEE::redDimPlotDefaults(data, 2)
  rdArgs$VisualBoxOpen <- c(TRUE, TRUE)
  rdArgs$ColorBy <- c("Column data", "Feature name")
  rdArgs$ColorByColData[1] <- "Annotation"
  
  feArgs <- iSEE::featAssayPlotDefaults(data, 1)
  
  feArgs$DataBoxOpen <- TRUE
  feArgs$XAxis <- "Column data"
  feArgs$XAxisColData <- "Annotation"
  
  color_fun <- function(n){
    c("gray","yellow","orange","red")
  }
  
  ecm <- iSEE::ExperimentColorMap(
    all_continuous=list( # shared
      assays=color_fun
    ),
    global_continuous=viridis::plasma # global
  )
  iSEE::iSEE(data, redDimArgs = rdArgs, featAssayArgs = feArgs, initialPanels=initialPanels, colormap = ecm)
}
