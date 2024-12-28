#BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                       'limma', 'S4Vectors', 'SingleCellExperiment',
#                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))

#### There is still connection problem when installing local github package.
#devtools::install_local('/BGFS1/projectdata/jasper/software/monocle3-master/')
####Timeout was reached: [api.github.com] Resolving timed out after 10000 milliseconds
####remove the last 3 rows in DESCRIPTION files:
#Remotes:
#    VPetukhov/ggrastr,
#    cole-trapnell-lab/leidenbase
####install depended github packages from local files:
#install.packages('ggrastr')
#devtools::install_local('/BGFS1/projectdata/jasper/software/leidenbase-master')
#### re-install monocle3 again. It worked.
##devtools::install_local('/BGFS1/projectdata/jasper/software/monocle3-master/')

#library(monocle3)
#message('Version of monocle3: ')
#message(packageVersion('monocle3'))

#' Run monocle 3
#' 
#' @param umap.min_dist Numeric indicating the minimum distance to be passed to
#'  UMAP function. Default is 0.1.See uwot package's umap for details. 
#'  umap.min_dist is smaller, points on UMAP are closer.
#'
#' @examples 
#' ## Don't run
#' run_monocle3(object = pbmc, annotation = 'Annotation',
#'  initial.celltype = 'NK',
#' use.seurat.preprocess = F, use.seurat.reduction = T,
#' seurat.preprocess.method = NULL, seurat.embeddings.method = 'umap',
#' save.files = T, outdir = './', return.obj = F, n.cores = 16)
#' mono3 <- readRDS('./mono3.rds')
#' plot_cells(mono3,color_cells_by = 'Annotation',
#' label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE)
#' 
#' @export
run_monocle3 <- function(object, 
                         annotation,
                         initial.celltype = NULL,
                         batch.effect = NULL,
                         batch.effect.continuous = NULL,
                         max_components = 2,
                         use.seurat.preprocess = F,
                         seurat.preprocess.method = NULL,
                         use.seurat.reduction = F,
                         seurat.embeddings.method = NULL,
                         num_dim = 50,
                         preprocess_method = 'PCA',
                         reduction_method = 'UMAP',
                         cluster_method = 'leiden',
                         k = 20,
                         umap.metric = "cosine",
                         umap.min_dist = 0.1,
                         umap.n_neighbors = 15L,
                         umap.fast_sgd = FALSE,
                         umap.nn_method = "annoy",
                         n.cores = 16,
                         use_partition = T,
                         close_loop = T,
                         learn_graph_control = NULL,
                         save.files = F, outdir = NULL,
                         return.obj = T) {
  if (!is.null(outdir)) {
    if (!dir.exists(outdir)) {
      message(outdir)
      stop('Not such directory.')
    }
  }
  expression_matrix <- object@assays$RNA@counts
  cell_metadata <- object@meta.data
  gene_annotation <- data.frame(id = rownames(object),
                                gene_short_name = rownames(object),
                                num_cells_expressed = apply(expression_matrix,1,
                                                            function(x) sum(x>0)))
  cds <- new_cell_data_set(expression_matrix,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_annotation)
  
  #Pre-process the data
  cds <- preprocess_cds(cds, num_dim = num_dim, method = preprocess_method)
  if (use.seurat.preprocess) {
    reducedDim(cds,type = preprocess_method) <- object@reductions[[seurat.preprocess.method]]@cell.embeddings
  }
  if (!is.null(batch.effect)) {
    cds <- align_cds(cds, 
                     alignment_group = batch.effect, 
                     residual_model_formula_str = batch.effect.continuous)
    preprocess_method <- 'Aligned'
  }
   #Reduce dimensionality and visualize the results
  cds <- reduce_dimension(cds,
                          max_components = max_components, 
                          preprocess_method = preprocess_method,
                          reduction_method = reduction_method,
                          umap.metric = umap.metric,
                          umap.min_dist = umap.min_dist,
                          umap.n_neighbors = umap.n_neighbors,
                          umap.fast_sgd = umap.fast_sgd,
                          umap.nn_method = umap.nn_method,
                          cores = n.cores
                          )
  if (use.seurat.reduction) {
    reducedDim(cds,type = reduction_method) <- object@reductions[[seurat.embeddings.method]]@cell.embeddings
  }
  #Cluster your cells
  if (use_partition) {
    cds <- cluster_cells(cds,
                         reduction_method = reduction_method, 
                         k = k, 
                         cluster_method = cluster_method)
    #plot_cells(cds, color_cells_by = "partition")
  }
  #Learn the trajectory graph
  cds <- learn_graph(cds, 
                     use_partition = use_partition, 
                     close_loop = close_loop, 
                     learn_graph_control = learn_graph_control)
  #cds <- order_cells(cds)
  # a helper function to identify the root principal points:
  if (is.null(initial.celltype)) {
    initial.celltype <- unique(colData(cds)[, annotation])[1]
  }
  get_earliest_principal_node <- function(cds, annotation, initial.celltype){
    cell_ids <- which(colData(cds)[, annotation] == initial.celltype)
    closest_vertex <-cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <- igraph::V(principal_graph(cds)[[reduction_method]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
    root_pr_nodes
  }
  cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds, annotation, initial.celltype))
  if (save.files) {
    saveRDS(cds, file = paste0(outdir,'mono3.rds'))
    p1 <- plot_cells(cds,
                    color_cells_by = annotation,
                    label_groups_by_cluster=FALSE,
                    label_leaves=FALSE,
                    label_branch_points=FALSE)
    p2 <- plot_cells(cds,
                    color_cells_by = "pseudotime",
                    label_cell_groups=FALSE,
                    label_leaves=FALSE,
                    label_branch_points=FALSE,
                    graph_label_size=1.5)
    pdf(paste0(outdir, 'monocle3_annotation.pdf'),8,6)
    print(p1)
    dev.off()
    pdf(paste0(outdir, 'monocle3_pseudotime.pdf'),8,6)
    print(p2)
    dev.off()
  }
  if (return.obj) {
    return(cds)
  }
}

