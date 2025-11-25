#' DEG using SCEPTRE
#' 
#' @param distance_threshold distance_threshold bases of one another and 
#' on the same chromosome are included in the cis discovery set.
#' @examples
#' object <- readRDS('../data/obj_pct01_logfc01.rds')
#' metaf <- object@@assays$RNA@@meta.features
#' grna_targetdf <- object@@assays$GDO@@meta.features
#' # select sgRNA id with position info
#' grna_targetdf <- grna_targetdf[(!is.na(grna_targetdf$chr)) | grepl('non-targeting',grna_targetdf$grna_target),]
#' grna_targetdf <- grna_targetdf[, c('grna_id', 'grna_target', 'chr', 'start', 'end')]
#' grna_targetdf <- left_join(grna_targetdf, data.frame(grna_target = metaf$response_name,
#'                                                    response_id = metaf$response_id))
#' grna_targetdf$grna_target[!is.na(grna_targetdf$response_id)] <- grna_targetdf$response_id[!is.na(grna_targetdf$response_id)]
#' grna_targetdf <- grna_targetdf[,-6]
#' grna_matrix <- object@@assays$GDO@@counts[grna_targetdf$grna_id,]
#' result_high <- DEG_SCEPTRE(hmoi, covariates = 'Sample', grna_targetdf, grna_matrix, moi = 'high')
#' result_low <- DEG_SCEPTRE(hmoi, covariates = 'Sample', grna_targetdf, grna_matrix, moi = 'low', grna_threshold = 2)
#' @export
DEG_SCEPTRE <- function(object, covariates, grna_targetdf, grna_matrix,
                        moi = "high",
                        to_ensembl = T,
                        control_group = 'default', 
                        p.adjust.method = 'fdr', 
                        grna_method = "default",
                        grna_threshold = 5,
                        distance_threshold = 5e6, p_mito_threshold = 0.075,
                        parallel = TRUE, n_processors = "auto") {
  if (to_ensembl) {
    # Transfering symbol to ensembl ID
    # response_matrix
    response_matrix <- object@assays$RNA@counts
    refgene <- gene_position_data_frame_grch38[!duplicated(gene_position_data_frame_grch38$response_name),]
    gene2ens <- left_join(data.frame(response_name = rownames(response_matrix)),
                          refgene)
    pos <- is.na(gene2ens$response_id)
    if (any(pos)) message('Genes that were not found in gene_position_data_frame_grch38: ', sum(pos))
    gene2ens <- gene2ens[!pos,]
    pos <- duplicated(gene2ens$response_id)
    if (any(pos)) message('Duplicated features were found when transfering symbol to ensembl ID. ', 'Remove genes: Ensembl ID:',  sum(pos))
    gene2ens <- gene2ens[!pos,]
    rownames(gene2ens) <- gene2ens$response_name
    response_matrix <- response_matrix[gene2ens$response_name,]
    rownames(response_matrix) <- gene2ens$response_id
  }
  
  print(dim(grna_matrix))
  #[1]   325 38856
  print(dim(response_matrix))
  #[1]   2468 38856
  print(head(grna_targetdf))
  sceptre_object <- import_data(
    response_matrix = response_matrix,
    grna_matrix = grna_matrix,
    grna_target_data_frame = grna_targetdf,
    moi = moi,
    extra_covariates = data.frame(batch = object@meta.data[, covariates]),
    response_names = gene2ens$response_name
  )
  positive_control_pairs <- construct_positive_control_pairs(
    sceptre_object = sceptre_object)
  
  if (moi == 'high') {
    discovery_pairs <- construct_cis_pairs(sceptre_object,
                                                   positive_control_pairs = positive_control_pairs,
                                                   distance_threshold = distance_threshold)
  }else{
    discovery_pairs <- construct_trans_pairs(
      sceptre_object = sceptre_object,
      positive_control_pairs = positive_control_pairs,
      pairs_to_exclude = "pc_pairs")
  }
 
  print(discovery_pairs)
  print(positive_control_pairs)

  #2-4. set analysis parameters, assign gRNAs, run qc, run calibration check
  
  sceptre_object <- sceptre_object |>
    set_analysis_parameters(
      discovery_pairs = discovery_pairs,
      positive_control_pairs = positive_control_pairs,
      side = "both",
      control_group = control_group,
      resampling_mechanism = "default",
      multiple_testing_method = p.adjust.method,
      multiple_testing_alpha = 0.999
  ) 
  print(sceptre_object)
  sceptre_object <- sceptre_object|>
    assign_grnas(parallel = parallel, n_processors = n_processors, 
                method = grna_method, threshold = grna_threshold) 
  
  sceptre_object <- sceptre_object |> 
    run_qc(
      p_mito_threshold = p_mito_threshold
    ) |>
    run_calibration_check(parallel = parallel,  n_processors = n_processors)
  
  
  # Running power check
  sceptre_object <- run_power_check(
    sceptre_object = sceptre_object,
    parallel = parallel, n_processors = n_processors
  )
  
  #  Run discovery analysis
  sceptre_object <- run_discovery_analysis(
    sceptre_object = sceptre_object,
    parallel = parallel, n_processors = n_processors
  )
  
  p <- plot(sceptre_object)
  
  discovery_result <- get_result(
    sceptre_object = sceptre_object,
    analysis = "run_discovery_analysis"
  )
  discovery_result <- left_join(discovery_result, metaf)
  return(list(object = sceptre_object, result = discovery_result, plot = p))
}
