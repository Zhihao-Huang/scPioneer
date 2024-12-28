#' Define outlier by MAD. Refer to https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html
#' 
#' @export
is_outlier <- function(metadata, metric, nmads = 5, mad.scale.factor = 1.0) {
  M <- metadata[,metric]
  coefMAD <- nmads * mad(M, constant = mad.scale.factor)
  outlier <- (M < median(M) - coefMAD) | (median(M) + coefMAD < M)
  return(outlier)
}

#' Blacklist is from doi.org/10.1038/s41586-022-05400-x
#' To avoid unexpected noise and expression artefacts by dissociation, a total of 1,514 genes associated with mitochondria (50 genes), heat-shock protein (178 genes), ribosome (1,253 genes) and dissociation (33 genes) were excluded (Supplementary Table 1).
#' 
#' @export
get_blacklist <- function(SeuratS4,param) {
  # change blastlist for regress if species is mouse. 
  blacklist <- blacklist
  if (param[['species']] == 'Mouse') {
    blacklistgene <- list()
    blacklistgene[['Mitochondria']] <- rownames(SeuratS4)[grepl("mt-", rownames(SeuratS4))]
    blacklistgene[['Ribosome']] <- rownames(SeuratS4)[grepl("^Rp[sl]", rownames(SeuratS4))]
    blacklistgene[['Heat.shock.protein']] <- rownames(SeuratS4)[grepl("^Hsp", rownames(SeuratS4))]
    df <- biospiper::convert_symbol_mouse_human(blacklist$Dissociation, 'human',
                                                using.local.file = T,
                                                file.dir = '/storage2/hlinglab/jasper/database/orthologs/')
    blacklistgene[['Dissociation']] <- df$MGI_symbol[!is.na(df$MGI_symbol)]
    maxl <- max(sapply(blacklistgene, length))
    n <- names(blacklistgene)
    blacklistgene <- lapply(blacklistgene, function(x) {
      if(length(x) < maxl) {
        c(x, rep(NA, maxl - length(x)))
      }else{x}
    })
    blacklist_mm10 <- do.call(cbind, blacklistgene)
    names(blacklist_mm10) <- n
    blacklist <- blacklist_mm10
  }
  return(blacklist)
}

#' Normalization module
#' 
#' @export
run_normalization <- function(SeuratS4, param) {
  blacklist <- get_blacklist(SeuratS4,param)
  noiselist <- as.list(blacklist)
  noiselist <- lapply(noiselist, function(x) x[x != ''])
  allgene <- as.vector(unlist(noiselist))
  noiselist[['noise_gene']] <- allgene
  geneset <- c("Mitochondria", "Heat.shock.protein", "Ribosome", "Dissociation","noise_gene")
  select_set <- 'noise_gene'
  regress_gene <- noiselist[select_set]
  # 5. Normalization
  if ('NULL' %in% param[['var.to.regress']]){
    vars.to.regress = NULL
  }else{
    vars.to.regress <- param[['var.to.regress']]
  }
  if (!'noise_gene' %in% colnames(SeuratS4@meta.data)) {
    # add noise gene score
    SeuratS4 <- AddModuleScore(SeuratS4, features = regress_gene,
                               name = 'noise_gene')
    colnames(SeuratS4@meta.data)[colnames(SeuratS4@meta.data) == 'noise_gene1'] <- 'noise_gene'
  }
  if (param[['normalize_method']] == 'LogNormalize') {
    assay <- 'RNA'
    SeuratS4 <- NormalizeData(object = SeuratS4, 
                              normalization.method = "LogNormalize",
                              scale.factor = param[["scale.factor"]])
    if (param[['remove_noise']] == 'TRUE') {
      SeuratS4 <- SeuratS4[!rownames(SeuratS4) %in% unlist(regress_gene),  ]
    }
    #############################################
    # 6. Detection of variable genes across the single cells
    SeuratS4 <- FindVariableFeatures(object = SeuratS4,
                                     selection.method = param[['vargene.method']],
                                     nfeatures = param[['nFeatures']])
    if (sum(param[['var.to.regress']] %in% select_set) > 0) {
      message('Excluding noise genes in high variable genes.')
      pos <- rownames(SeuratS4) %in% unlist(regress_gene)
      SeuratS4@assays$RNA@meta.data$var.features[pos] <- NA
    }
    if (sum(param[['var.to.regress']] %in% c('S.Score','G2M.Score')) > 1) {
      message('Excluding cellcycle genes in high variable genes.')
      pos <- rownames(SeuratS4) %in% unlist(cc.genes)
      SeuratS4@assays$RNA@meta.data$var.features[pos] <- NA
    }
    # Identify the 10 most highly variable genes
    top10 <- head(VariableFeatures(SeuratS4), 10)
    # plot variable features with and without labels
    plot1 <- VariableFeaturePlot(SeuratS4)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    pdf(paste(param[['outdir']],"Filter_genes.pdf",sep=''),width = 8,height = 7)
    print(plot1)
    print(plot2)
    dev.off()
    
    ############################################
    # 7. Scaling the data
    message('7. Scaling the data.')
    ###regress
    if (param[['scale.all.gene']]) {
      scale.genes <- rownames(SeuratS4)
    }else{
      scale.genes <- VariableFeatures(SeuratS4)
    }
    SeuratS4 <- ScaleData(SeuratS4,features = scale.genes, 
                          vars.to.regress = vars.to.regress)
  }else if (param[['normalize_method']] == 'SCT') {
    assay <- 'SCT'
    if (param[['remove_noise']] == 'TRUE') {
      SeuratS4 <- SeuratS4[!rownames(SeuratS4) %in% unlist(regress_gene),  ]
    }
    SeuratS4 <- SCTransform(SeuratS4, vars.to.regress = vars.to.regress, 
                            variable.features.n = param[['nFeatures']])
    if (param[['var.to.regress']] %in% select_set) {
      pos <- rownames(SeuratS4) %in% unlist(regress_gene)
      SeuratS4@assays$SCT@var.features <- SeuratS4@assays$SCT@var.features[!pos] 
    }
    if (sum(param[['var.to.regress']] %in% c('S.Score','G2M.Score')) > 1) {
      message('Excluding cellcycle genes in high variable genes.')
      pos <- rownames(SeuratS4) %in% unlist(regress_gene)
      SeuratS4@assays$SCT@var.features <- SeuratS4@assays$SCT@var.features[!pos] 
    }
  }else{
    stop('Only support normalization methods: LogNormalize or SCT.')
  }
  return(SeuratS4)
}

#' Running PCA
#' 
#' @export
run_pca <- function(SeuratS4, param) {
  SeuratS4 <- RunPCA(SeuratS4, npcs = param[['npcs']], 
                     features = VariableFeatures(object = SeuratS4))
  p <- ElbowPlot(SeuratS4, ndims = param[['npcs']])
  ggsave(paste(param[['outdir']],"Select_pc.pdf",sep=''), p, width = 8,height = 7)
  p <- DimPlot(SeuratS4, reduction = "pca",group.by = 'Sample')
  ggsave(paste(param[['outdir']],"pca.pdf",sep=''), p, width = 8,height = 7)
  return(SeuratS4)
}

#' Detect doublet and filter low-quality cells
#' 
#' @export
run_doublet_and_filter <- function(SeuratS4, param) {
  if (param[['detect.doublet']] == 'scrublet') {
    ###########################################
    # scrublet
    message('3. Detect doublet by scrublet.')
    # get raw meta data for statistic
    meta <- SeuratS4@meta.data[,c("Sample","percent.mt", "percent.ribo",
                                  "nFeature_RNA","nCount_RNA")]
    # filter and normalize for detecting doublet.
    if (!'NULL' %in% param[['mad.outlier.metric']]) {
      metalist <- list()
      for (sam in unique(SeuratS4$Sample)) {
        subobj <- subset(SeuratS4, subset = Sample == sam)
        #subobj <- subset(SeuratS4, Sample == sam)
        metadata <- subobj@meta.data
        poslist <- list()
        metricdf <- data.frame(param[['mad.outlier.metric']], param[['mad.outlier.coef']], param[['mad.outlier.constant']])
        for (i in 1:nrow(metricdf)) poslist[[i]] <- is_outlier(metadata, metricdf[i,1], nmads = metricdf[i,2], mad.scale.factor = metricdf[i,3])
        metadata$mad.outliers <- Reduce('|', poslist)
        message(paste0('Remove MAD outlier of Sample ',sam,': ',sum(metadata$mad.outliers)))
        metalist[[sam]] <- metadata
      }
      metadata <- do.call(rbind, metalist)
      rownames(metadata) <- metadata$Cellname
      metadata <- metadata[colnames(SeuratS4),]
      SeuratS4$mad.outliers <- metadata$mad.outliers
      SeuratS4 <- subset(SeuratS4, subset = mad.outliers == FALSE)
    }
    SeuratS4 <- subset(SeuratS4, subset = nFeature_RNA > param[['min.features']] & 
                         nFeature_RNA < param[['max.features']] & 
                         nCount_RNA > param[['min.UMIs']] &
                         nCount_RNA < param[['max.UMIs']] &
                         percent.mt < param[['max.mt']])
    
    SeuratS4 <- run_normalization(SeuratS4, param)
    SeuratS4 <- RunPCA(SeuratS4)
    SeuratS4 <- RunUMAP(SeuratS4, reduction = param[['reduction_name']],
                        reduction.name = param[['umap_name']],
                        dims =  param[['select_PCs']])
    # Seurat to CellDataSet object
    #message('Convert Seurat object to monocle2 object...')
    #cds <- as.CellDataSet(SeuratS4, assay = 'RNA', reduction = 'umap')
    message("Running scrublet...")
    scrublet_res <- scrublet_R(SeuratS4, python_home = param[['python_home']], 
                               return_results_only = T,
                               scrublet_script_path = param[['scrublet_script_path']],
                               min_cells = param[['min.cells']],
                               expected_doublet_rate = param[['expected.doublet.rate']])
    message("Scrublet done.")
    SeuratS4$doublet_scores <- scrublet_res$doublet_scores
    p <- ggplot(data.frame(doublet_scores = SeuratS4$doublet_scores), 
                aes(x = doublet_scores)) + geom_histogram(bins = 200) + 
      scale_y_log10() + theme_bw()
    p <- p + geom_vline(xintercept = param[['max.doublet.rate']],
                        linetype="dashed", color = "black", size=1)
    ggsave(paste0(param[['outdir']], 'hist_doublet.pdf'),p, width = 6, height = 3.2)
    dbpos <- SeuratS4$doublet_scores > param[['max.doublet.rate']]
    sum(dbpos)
    #[1] 1529
    SeuratS4$doublet <- 'singlet'
    SeuratS4$doublet[dbpos] <- 'doublet'
    p1 <- FeaturePlot2(SeuratS4, features = 'doublet_scores')
    p2 <- DimPlot(SeuratS4, group.by = 'doublet')
    SeuratMajorVersion <- gsub('\\..*$','',packageVersion('Seurat'))
    device <- ifelse(SeuratMajorVersion == 4, yes = 'pdf', no = 'png')
    ggsave(paste0(param[['outdir']],'doublet_score.',device), p1 + p2, device = device,
           width = 7.7, height = 3.9)
    saveRDS(SeuratS4, file = paste0(param[['outdir']],'rawObject.rds'))
    ##raw statistic. add doublet info. to raw meta data
    meta$Cellname <- rownames(meta)
    meta <- left_join(meta, SeuratS4@meta.data[,c("Cellname","doublet")], by = 'Cellname')
    rownames(meta) <- meta$Cellname
    meta$doublet[is.na(meta$doublet)] <- 'doublet'
    meta$Cellname <- NULL
    statmat <- stat_num(meta, min_nFeature_RNA = param[['min.features']],
                        max_nFeature_RNA = param[['max.features']],
                        percent.mt = param[['max.mt']],
                        str.param.list = list(doublet = 'singlet'))
    write.table(statmat, file = paste0(param[['outdir']], 'sample_cell_number.txt'), 
                quote = F, sep = '\t')
    print(dim(SeuratS4))
    #[1] 23467 56089
    if (param[['filter.doublet']]) {
      #SeuratS4 <- readRDS(paste0(outdir,'rawObject.rds'))
      SeuratS4 <- subset(SeuratS4, subset = nFeature_RNA > param[['min.features']] & 
                           nFeature_RNA < param[['max.features']] &
                           percent.mt < param[['max.mt']] & 
                           doublet_scores <= param[['max.doublet.rate']]) 
    }
    
  }else if (param[['detect.doublet']] == 'scDblFinder') {
    message('3. Detect doublet by scDblFinder.')
    # get raw meta data for statistic
    meta <- SeuratS4@meta.data[,c("Sample","percent.mt", "percent.ribo",
                                  "nFeature_RNA","nCount_RNA")]
    # get clusters
    if (!'NULL' %in% param[['mad.outlier.metric']]) {
      metalist <- list()
      for (sam in unique(SeuratS4$Sample)) {
        subobj <- subset(SeuratS4, subset = Sample == sam)
        #subobj <- subset(SeuratS4, Sample == sam)
        metadata <- subobj@meta.data
        poslist <- list()
        metricdf <- data.frame(param[['mad.outlier.metric']], param[['mad.outlier.coef']], param[['mad.outlier.constant']])
        for (i in 1:nrow(metricdf)) poslist[[i]] <- is_outlier(metadata, metricdf[i,1], nmads = metricdf[i,2], mad.scale.factor = metricdf[i,3])
        metadata$mad.outliers <- Reduce('|', poslist)
        message(paste0('Remove MAD outlier of Sample ',sam,': ',sum(metadata$mad.outliers)))
        metalist[[sam]] <- metadata
      }
      metadata <- do.call(rbind, metalist)
      rownames(metadata) <- metadata$Cellname
      metadata <- metadata[colnames(SeuratS4),]
      SeuratS4$mad.outliers <- metadata$mad.outliers
      SeuratS4 <- subset(SeuratS4, subset = mad.outliers == FALSE)
    }
    SeuratS4 <- subset(SeuratS4, subset = nFeature_RNA > param[['min.features']] & 
                         nFeature_RNA < param[['max.features']] & 
                         nCount_RNA > param[['min.UMIs']] &
                         nCount_RNA < param[['max.UMIs']] &
                         percent.mt < param[['max.mt']])
    
    SeuratS4 <- run_normalization(SeuratS4, param)
    SeuratS4 <- RunPCA(SeuratS4)
    set.seed(123)
    SeuratS4 <- RunUMAP(SeuratS4, dims =  param[['select_PCs']],
                        reduction = param[['reduction_name']],
                        reduction.name = param[['umap_name']])
    SeuratS4 <- FindNeighbors(SeuratS4, dims = param[['select_PCs']],
                              reduction = param[['reduction_name']], 
                              k.param = param[['k.param']])
    # resolution was set to low value for remaining main clusters
    SeuratS4 <- FindClusters(SeuratS4, resolution = 0.3,
                             algorithm = param[['cluster_algorithm']])
    # detect doublet using scDblFinder
    message("Running scDblFinder...")
    if ( param[['scdblFinder_dbr']] == 'NULL') param[['scdblFinder_dbr']] <- NULL
    if ( param[['scdblFinder_dbr.sd']] == 'NULL') param[['scdblFinder_dbr.sd']] <- NULL
    if (param[['is_multidata']] == 'TRUE') {
      bp <- BiocParallel::MulticoreParam(4, RNGseed=123)
      sce <- scDblFinder::scDblFinder(GetAssayData(SeuratS4, layer="counts"), 
                                      clusters = Idents(SeuratS4),
                                      samples = SeuratS4$Sample,
                                      BPPARAM = bp,
                                      dbr = param[['scdblFinder_dbr']],
                                      dbr.sd = param[['scdblFinder_dbr.sd']]
      )
    }else{
      set.seed(123)
      sce <- scDblFinder::scDblFinder(GetAssayData(SeuratS4, layer="counts"), 
                                      clusters = Idents(SeuratS4),
                                      dbr = param[['scdblFinder_dbr']],
                                      dbr.sd = param[['scdblFinder_dbr.sd']]
      )
    }
    message("scDblFinder running is over.")
    # port the resulting scores back to the Seurat object:
    SeuratS4$doublet_scores <- sce$scDblFinder.score
    SeuratS4$doublet <- sce$scDblFinder.class
    p1 <- FeaturePlot2(SeuratS4, features = 'doublet_scores')
    p2 <- DimPlot(SeuratS4, group.by = 'doublet')
    SeuratMajorVersion <- gsub('\\..*$','',packageVersion('Seurat'))
    device <- ifelse(SeuratMajorVersion == 4, yes = 'pdf', no = 'png')
    ggsave(paste0(param[['outdir']],'doublet_score.',device), p1 + p2, device = device,
           width = 7.7, height = 3.9)
    saveRDS(SeuratS4, file = paste0(param[['outdir']],'rawObject.rds'))
    ##raw statistic. add doublet info. to raw meta data
    meta$Cellname <- rownames(meta)
    SeuratS4$Cellname <- rownames(SeuratS4@meta.data)
    meta <- left_join(meta, SeuratS4@meta.data[,c("Cellname","doublet")], by = 'Cellname')
    rownames(meta) <- meta$Cellname
    meta$doublet[is.na(meta$doublet)] <- 'doublet'
    meta$Cellname <- NULL
    statmat <- stat_num(meta, min_nFeature_RNA = param[['min.features']],
                        max_nFeature_RNA = param[['max.features']],
                        percent.mt = param[['max.mt']],
                        str.param.list = list(doublet = 'singlet'))
    write.table(statmat, file = paste0(param[['outdir']], 'sample_cell_number.txt'), 
                quote = F, sep = '\t')
    print(dim(SeuratS4))
    if (param[['filter.doublet']]) {
      #SeuratS4 <- readRDS(paste0(outdir,'rawObject.rds'))
      SeuratS4 <- subset(SeuratS4, subset = nFeature_RNA > param[['min.features']] & 
                           nFeature_RNA < param[['max.features']] &
                           percent.mt < param[['max.mt']] & 
                           doublet == 'singlet') 
    }
  }
  # not detect doublet
  else{
    ##raw statistic
    meta <- SeuratS4@meta.data[,c("Sample","percent.mt", "percent.ribo",
                                  "nFeature_RNA","nCount_RNA")]
    statmat <- stat_num(meta, min_nFeature_RNA = param[['min.features']],
                        max_nFeature_RNA = param[['max.features']],
                        percent.mt = param[['max.mt']])
    write.table(statmat, file = paste0(param[['outdir']], 'sample_cell_number.txt'), 
                quote = F, sep = '\t')
    print(dim(SeuratS4))
    if (!'NULL' %in% param[['mad.outlier.metric']]) {
      metalist <- list()
      for (sam in unique(SeuratS4$Sample)) {
        subobj <- subset(SeuratS4, subset = Sample == sam)
        #subobj <- subset(SeuratS4, Sample == sam)
        metadata <- subobj@meta.data
        poslist <- list()
        metricdf <- data.frame(param[['mad.outlier.metric']], param[['mad.outlier.coef']], param[['mad.outlier.constant']])
        for (i in 1:nrow(metricdf)) poslist[[i]] <- is_outlier(metadata, metricdf[i,1], nmads = metricdf[i,2], mad.scale.factor = metricdf[i,3])
        metadata$mad.outliers <- Reduce('|', poslist)
        message(paste0('Remove MAD outlier of Sample ',sam,': ',sum(metadata$mad.outliers)))
        metalist[[sam]] <- metadata
      }
      metadata <- do.call(rbind, metalist)
      rownames(metadata) <- metadata$Cellname
      metadata <- metadata[colnames(SeuratS4),]
      SeuratS4$mad.outliers <- metadata$mad.outliers
      SeuratS4 <- subset(SeuratS4, subset = mad.outliers == FALSE)
    }
    SeuratS4 <- subset(SeuratS4, subset = nFeature_RNA > param[['min.features']] & 
                         nFeature_RNA < param[['max.features']] & 
                         nCount_RNA > param[['min.UMIs']] &
                         nCount_RNA < param[['max.UMIs']] &
                         percent.mt < param[['max.mt']])
  }
  return(SeuratS4)
}

#' Run Batch effect correction based on Incorporated methods in Seurat V5
#' 
#' @export
run_batch_effect_correction <- function(SeuratS4, param) {
  if (param[['batch_correct_method']] != 'NULL') {
    if (param[['normalize_method']] == 'LogNormalize') {
      assay = 'RNA'
    }else if (param[['normalize_method']] == 'SCT') {
      assay = 'SCT'
    }
    SeuratS4[[assay]] <- split(SeuratS4[[assay]], f = SeuratS4@meta.data[,param[['batch_correct_colname']]])
    SeuratS4 <- run_normalization(SeuratS4, param)
    SeuratS4 <- RunPCA(SeuratS4)
    SeuratS4 <- FindNeighbors(SeuratS4, dims = param$select_PCs, reduction = "pca")
    SeuratS4 <- FindClusters(SeuratS4, resolution = param$resolution, cluster.name = "unintegrated_clusters")
    SeuratS4 <- RunUMAP(SeuratS4, dims = param$select_PCs, reduction = "pca", reduction.name = "umap.unintegrated")
    if (param[['sketch_integration']] %in% c('TRUE','T')) {
      # Perform integration on the sketched cells across samples
      SeuratS4 <- SketchData(object = SeuratS4, ncells = param[['sketch_ncells']], 
                             method = "LeverageScore", sketched.assay = "sketch")
      DefaultAssay(SeuratS4) <- "sketch"
      SeuratS4 <- run_normalization(SeuratS4, param)
      SeuratS4 <- RunPCA(SeuratS4)
    } 
    # SeuratV3-CCA and harmony are available
    if (param[['batch_correct_method']] == 'cca') {
      message('Perform integraton by CCA')
      # a. SeuratV5-CCA
      SeuratS4 <- IntegrateLayers(
        object = SeuratS4, method = CCAIntegration,
        orig.reduction = "pca", new.reduction = "integrated.cca",
        verbose = FALSE, dims = param[['select_PCs']],
        normalization.method = param[['normalize_method']],
        k.weight = param[['k.weight']]
      )
      param[['reduction_name']] <- 'integrated.cca'
      param[['umap_name']] <- 'umap.CCA'
    }
    if (param[['batch_correct_method']] == 'harmony') {
      message('Perform integraton by harmony')
      # b. harmony
      SeuratS4 <- IntegrateLayers(
        object = SeuratS4, method = HarmonyIntegration,
        orig.reduction = "pca", new.reduction = "harmony",
        verbose = FALSE, dims = param[['select_PCs']],
        normalization.method = param[['normalize_method']]
      )
      param[['reduction_name']] <- 'harmony'
      param[['umap_name']] <- 'umap.Harmony'
    }
    if (param[['batch_correct_method']] == 'rpca') {
      message('Perform integraton by rpca')
      # b. RPCA
      SeuratS4 <- IntegrateLayers(
        object = SeuratS4, method = RPCAIntegration,
        orig.reduction = "pca", new.reduction = "integrated.rpca",
        verbose = FALSE, dims = param[['select_PCs']],
        normalization.method = param[['normalize_method']],
        k.weight = param[['k.weight']]
      )
      param[['reduction_name']] <- 'integrated.rpca'
      param[['umap_name']] <- 'umap.RPCA'
    }
    if (param[['batch_correct_method']] == 'mnn') {
      message('Perform integraton by mnn')
      # c. FastMNN
      SeuratS4 <- IntegrateLayers(
        object = SeuratS4, method = FastMNNIntegration,
        orig.reduction = "pca", new.reduction = "integrated.mnn",
        verbose = FALSE
      )
      param[['reduction_name']] <- 'integrated.mnn'
      param[['umap_name']] <- 'umap.MNN'
    }
    if (param[['batch_correct_method']] == 'scvi') {
      message('Perform integraton by scvi')
      # Must use specific version of Seurat before running scvi. See https://github.com/satijalab/seurat/issues/7944
      SeuratS4 <- IntegrateLayers(
        object = SeuratS4, method = scVIIntegration, 
        groups = "Method", dims = param[['select_PCs']],
        normalization.method = param[['normalize_method']],
        orig.reduction = "pca", new.reduction = "integrated.scvi",
        conda_env = param[['python_env']], verbose = F
      )
      param[['reduction_name']] <- 'integrated.scvi'
      param[['umap_name']] <- 'umap.SCVI'
    }
    if (param[['sketch_integration']] %in% c('TRUE','T')) {
      # Integrate the full datasets based on sketch cells.The non-sketched cells are not loaded into memory.
      object <- ProjectIntegration(object = object, sketched.assay = "sketch", 
                                   assay = assay, reduction = param[['reduction_name']])
      param[['reduction_name']] <- paste0(param[['reduction_name']],".full")
      object <- ProjectData(object = object, sketched.assay = "sketch", assay = assay,
                            sketched.reduction = param[['reduction_name']],
                            full.reduction = param[['reduction_name']],
                            dims = param[['select_PCs']])
      param[['umap_name']] <- paste0(param[['umap_name']],'.Full')
    }
    SeuratS4 <- JoinLayers(SeuratS4)
  }else{
    message('Skip batch effect correction.')
  }
  return(list(obj = SeuratS4, par = param))
}


#' Run QC and clustering using Seurat v5
#' 
#' @export
PHASE1_run_Seurat_v5_QC_clustering <- function(param) {
  ##########################################
  ##### check parameter list
  
  if (class(param) != 'list') {
    stop('Please input parameter using list object.')
  }
  template <- PHASE1_run_Seurat_v5_QC_clustering_param_template()
  missing <- names(template)[!names(template) %in% names(param)]
  if (length(missing) > 0) {
    message('Some parameters were not specified, and would be set to default as showed below.')
    print(template[names(template) %in% missing])
  }
  param[missing] <- template[missing]
  #message('All parameters:')
  #print(param)
  ###########################################
  if (!dir.exists(param[['outdir']])) dir.create(param[['outdir']])
  if (!dir.exists(param[['outdir']])) stop('Could not create out directory: ',param[['outdir']])
  ##########################################
  
  
  ##########################################
  # 1. Create Seurat object
  message('#########################  1. Create Seurat object.  #########################')
  if (param[['object']] == 'NULL') {
    message('Load raw metrics from:')
    message(param[['samplelist']])
    samplelist <- read.table(param[['samplelist']], sep = '\t',header = T, stringsAsFactors = F)
    print(head(samplelist))
    objlist <- list()
    for (i in 1:nrow(samplelist)) {
      data <- Read10X(data.dir = samplelist$datadir[i])
      dataset_name = samplelist$samplename[i]
      colnames(data) <- paste0(dataset_name,'_cell_',colnames(data))
      obj <- CreateSeuratObject(counts=data, min.cells = 0, 
                                min.features=0, project=param[['project']])
      objlist[[dataset_name]] <- obj
    }
    SeuratS4 <- Reduce(merge, objlist)
    SeuratS4$Cellname <- colnames(SeuratS4)
    if (param[['is_multidata']] == 'TRUE') {
      SeuratS4@meta.data$Sample <- gsub('_cell_.*$','',colnames(SeuratS4))
    }else{
      SeuratS4@meta.data$Sample <- param[['project']]
    }
  }else{
    # single sample
    main_obj <- readRDS(param[['object']])
    if (param[['subset_colname']] != 'NULL') {
      cellID <- colnames(main_obj)[main_obj@meta.data[,param[['subset_colname']]] %in% 
                                     param[['subset_type']]]
      SeuratS4 <- subset(main_obj, cells = cellID)
      message(paste0('subset ',paste(param[['subset_type']],collapse = ','),'(',length(cellID), ')',
                     ' from main object (', ncol(main_obj),')'))
    }else{
      SeuratS4 <- main_obj
    }
    SeuratS4 <- CreateSeuratObject(counts = GetAssayData(SeuratS4, assay = 'RNA', layer = 'counts'), 
                                   meta.data = SeuratS4@meta.data,
                                   min.cells = param[['min.cells']],
                                   project=param[['project']])
  }
  if (!'Sample' %in% colnames(SeuratS4@meta.data)) {
    if (param[['sample_colname']] != 'NULL') {
      SeuratS4@meta.data$Sample <- SeuratS4@meta.data[,param[['sample_colname']]]
    }else{
      message('Sample was not specified in meta.data. Use org.ident instead.')
      SeuratS4@meta.data$Sample <- SeuratS4@meta.data$orig.ident
    }
  }
  SeuratS4$Sample <- as.vector(SeuratS4$Sample)
  SeuratS4$Cellname <- colnames(SeuratS4)
  SeuratS4$log1pCount <- log(SeuratS4$nCount_RNA + 1)
  SeuratS4$log1pGene <- log(SeuratS4$nFeature_RNA + 1)
  ############################################
  # 2. QC of raw data
  message('#########################  2.Run QC before filtering.  #########################')
  #SeuratS4 <- QC_raw(SeuratS4,outdir,species,multi.samples=is_multidata)
  message(param[['species']])
  SeuratS4 <- QC_raw(object = SeuratS4, species = param[['species']],
                     outdir = param[['outdir']], 
                     max.mt = param[['max.mt']],
                     min.features = param[['min.features']], 
                     max.features = param[['max.features']],
                     multi.samples = param[['is_multidata']], 
                     do_cellcycle = param[['cal_cellcycle']],
                     plot.raw = param[['plot.raw']])
  
  
  ###########################################
  # 3. Detect doublet 
  message('#########################  3. Find low-quality cells and filter.  #########################')
  SeuratS4 <- run_doublet_and_filter(SeuratS4, param)
 
  #4. QC plot after filtering
  message('#########################  4. Run QC after filtering.  #########################')
  message('Complete filtering.')
  print(dim(SeuratS4))
  num_g <- length(unique(SeuratS4@meta.data$Sample))
  features1 = c('percent.mt','nFeature_RNA',"nCount_RNA",
                "percent.ribo", "percent.HSP", 
                "S.Score","G2M.Score")
  #p <- RidgePlot(SeuratS4, features = features1, ncol = 1,group.by = 'Sample') + ylab('')+
  #  theme(legend.position = 'none')
  # ggsave(paste(outdir, "QC.RidgePlot_density_filtered.pdf", sep = ""),p,
  #        width = 6, height = (8 + num_g * 1.2))
  
  #################################################
  # 5. Normalization
  message('#########################  5. Perform normalization. #########################')
  SeuratS4 <- run_normalization(SeuratS4, param)
  
  ###########################################
  # 8. Perform linear dimensional reduction (PCA)
  message('#########################  6. Perform linear dimensional reduction (PCA). ########################')
  SeuratS4 <- run_pca(SeuratS4, param)
  
  ##############################################
  # 9. Batch effect correction
  message('#########################  7. Batch effect correction.  #########################')
  obj_param <- run_batch_effect_correction(SeuratS4, param)
  SeuratS4 <- obj_param$obj
  param <- obj_param$par
  rm(obj_param)
  ###############################################
  # 10. Run Non-linear dimensional reduction
  message('#########################  8. Run Non-linear dimensional reduction.  #########################')
  set.seed(123)
  SeuratS4 <- RunUMAP(SeuratS4, reduction = param[['reduction_name']],
                      reduction.name = param[['umap_name']],
                      dims =  param[['select_PCs']])
  SeuratS4@meta.data$UMAP_1 <- SeuratS4@reductions[[param[['umap_name']]]]@cell.embeddings[,1]
  SeuratS4@meta.data$UMAP_2 <- SeuratS4@reductions[[param[['umap_name']]]]@cell.embeddings[,2]
  
  #############################################
  # 11. Cluster the cells
  message('#########################  9. Cluster the cells.  #########################')
  set.seed(123)
  SeuratS4 <- FindNeighbors(SeuratS4, dims = param[['select_PCs']],
                            reduction = param[['reduction_name']], 
                            k.param = param[['k.param']])
  SeuratS4 <- FindClusters(SeuratS4, resolution = param[['resolution']],
                           algorithm = param[['cluster_algorithm']])
  p1 <- DimPlot(SeuratS4, reduction = param[['umap_name']])
  p2 <- DimPlot(SeuratS4, group.by = 'Sample', reduction = param[['umap_name']])
  ggsave(paste0(param[['outdir']], 'UMAP.pdf'), patchwork::wrap_plots(p1,p2), width = 12, height = 4.8)
  ############################################
  # 12. Do DEG analysis of raw clusters 
  if (param[['do_DEG_wilcox']] == 'TRUE') {
    message('#########################  10. Do DEG analysis of raw clusters.  #########################')
    mmk <- FindAllMarkers(SeuratS4, only.pos = T)
    write.table(mmk, quote = F, sep = '\t', 
                file = paste0(param[['outdir']], 'FindAllMarkers.xls'))
  }
  ############################################
  return(SeuratS4)
}

#' Parameter template for PHASE1_run_Seurat_v3_QC_clustering.
#' 
#' @export
PHASE1_run_Seurat_v5_QC_clustering_param_template <- function() {
  ##config
  param <- list()
  param[['samplelist']] <- "NULL"
  param[['object']] <- "./result/101.cluster_main/final.rds"
  param[['subset_colname']] <- 'NULL'
  param[['subset_type']] <- c('Macrophage')
  param[['sample_colname']] <- 'NULL'
  param[['project']] <- 'Example_project'
  param[['outdir']] <- './result/'
  param[['min.cells']] <- 0
  param[['min.features']] <- 200
  param[['max.features']] <- Inf
  param[['min.UMIs']]  <- -Inf
  param[['max.UMIs']] <- Inf
  param[['max.mt']] <- 10
  param[['mad.outlier.metric']] <- 'NULL'
  param[['mad.outlier.coef']] <- 5
  param[['mad.outlier.constant']] <- 1
  param[['cal_cellcycle']] <- 'TRUE'
  param[['plot.raw']] <- 'TRUE'
  param[['normalize_method']] <- 'LogNormalize' #### 'LogNormalize','SCT'
  param[["scale.factor"]] = 10000
  param[['var.to.regress']] <- c('NULL')
  param[['scale.all.gene']] <- FALSE
  param[['nFeatures']] <- 2000
  param[['vargene.method']] = 'vst'
  param[['tsne_method']] <- 'Rtsne'
  param[['npcs']] = 50
  param[['dims.use']] <- 30
  param[['select_PCs']] <- 1:param[['dims.use']]
  param[['reduction_name']] <- 'pca'
  param[['batch_correct_method']] <- 'NULL' #### 'cca','rpca','harmony','fastmnn','scvi'
  param[['batch_correct_colname']] <- 'Sample'
  param[['k.param']] <- 30
  param[['resolution']] <- 0.8
  param[['cluster_algorithm']] <- 1
  param[['perplexity']] <- 30
  param[['nn.eps']] <- 0
  param[['n.start']] <- 50
  param[['min_mean']] <- 0.0125
  param[['max_mean']] <- 3
  param[['min_disp']] <- 0.5
  param[['detect.doublet']] <- 'NULL' ### scDblFinder or scrublet
  param[['scdblFinder_dbr']] <- 'NULL'
  param[['scdblFinder_dbr.sd']] <- 'NULL'
  param[['filter.doublet']] <- F
  param[['python_home']] <- '/home/jasper/.conda/envs/Seurat_v5/bin/python'
  param[['python_env']] <- gsub('bin.*$','',param[['python_home']])
  param[['scrublet_script_path']] <- '/storage2/hlinglab/jasper/pipeline_development/scPioneer/scpioneer_v2.0.0/scPioneer/R/scrublet.py'
  param[['expected.doublet.rate']] <- 0.05
  param[['max.doublet.rate']] <- 0.2
  param[['is_multidata']] <- 'TRUE'
  param[['species']] <- 'Human'
  param[['do_DEG_wilcox']] <- 'FALSE'
  param[['remove_noise']] <- 'FALSE'
  param[['k.weight']] <- 100
  param[['umap_name']] <- 'umap'
  param[['sketch_integration']] <- 'FALSE'
  param[['sketch_ncells']] <- 5000
  return(param)
}

