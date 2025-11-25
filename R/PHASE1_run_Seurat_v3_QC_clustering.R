#' Run QC and clustering using Seurat v3.
#' 
#' @export
PHASE1_run_Seurat_v3_QC_clustering <- function(param) {
  ##########################################
  ##### check parameter list
  
  if (class(param) != 'list') {
    stop('Please input parameter using list object.')
  }
  template <- PHASE1_run_Seurat_v3_QC_clustering_param_template()
  missing <- names(template)[!names(template) %in% names(param)]
  if (length(missing) > 0) {
    message('Some parameters were not specified, and would be set to default as showed below.')
    print(template[names(template) %in% missing])
  }
  param[missing] <- template[missing]
  #message('All parameters:')
  #print(param)
  ###########################################
  outdir <- param[['outdir']]
  if (!dir.exists(outdir)) dir.create(outdir)
  if (!dir.exists(outdir)) stop('Could not create out directory: ',outdir)
  ##########################################
  
  
  ##########################################
  # 1. Create Seurat object
  message('1. Create Seurat object.')
  if (param[['object']] == 'NULL') {
    message('Load raw metrics from:')
    message(param[['samplelist']])
    samplelist <- read.table(param[['samplelist']], sep = '\t',header = T, stringsAsFactors = F)
    print(head(samplelist))
    # Two methods for merging raw data
    # 1. merge raw data then create object
    oldmethod <- F
    if (oldmethod) {
      datalist <- list()
      for (i in 1:nrow(samplelist)) {
        data <- Read10X(data.dir = samplelist$datadir[i])
        dataset_name = samplelist$samplename[i]
        colnames(data) <- paste0(dataset_name,'_cell_',colnames(data))
        datalist[[dataset_name]] <- data
      }
      data <- do.call(cbind,datalist)
      #Convert '-' to '.' in cell names.
      colnames(data) <- gsub('-','\\.',colnames(data))
      SeuratS4 <- CreateSeuratObject( counts=data, min.cells = param[['min.cells']], 
                                      min.features=0, project=param[['project']])
    }
    # 2. create objects then merge objects
    mergeObj <- T
    if (mergeObj) {
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
    }
    
    SeuratS4$Cellname <- colnames(SeuratS4)
    if (param[['is_multidata']] == 'TRUE') {
      SeuratS4@meta.data$Sample <- gsub('_cell_.*$','',colnames(SeuratS4))
    }else{
      SeuratS4@meta.data$Sample <- param[['project']]
    }
    
  }else{
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
    SeuratS4 <- CreateSeuratObject(counts = SeuratS4@assays$RNA@counts, 
                                   meta.data = SeuratS4@meta.data,
                                   min.cells = param[['min.cells']],
                                   project=param[['project']])
  }
  if (param[['sample_colname']] != 'NULL') {
    SeuratS4@meta.data$Sample <- SeuratS4@meta.data[,param[['sample_colname']]]
  }
  SeuratS4$Sample <- as.vector(SeuratS4$Sample)
  ############################################
  # 2. QC of raw data
  message('2.Running QC before filtering.')
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
  # 4. Filter cells.
  message('4. Filter cells.')
  # change blastlist for regress if species is mouse. 
  blacklist <- blacklist
  if (param[['species']] == 'Mouse') {
    blacklistgene <- list()
    blacklistgene[['Mitochondria']] <- rownames(SeuratS4)[grepl("mt-", rownames(SeuratS4))]
    blacklistgene[['Ribosome']] <- rownames(SeuratS4)[grepl("^Rp[sl]", rownames(SeuratS4))]
    blacklistgene[['Heat.shock.protein']] <- rownames(SeuratS4)[grepl("^Hsp", rownames(SeuratS4))]
    df <- biospiper::convert_symbol_mouse_human(blacklist$Dissociation, 'human',using.local.file = T,file.dir = '/storage2/hlinglab/jasper/database/orthologs/')
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
  if (param[['detect.doublet']] == 'scrublet') {
    ###########################################
    # 3. Detect doublet by scrublet
    message('3. Detect doublet by scrublet.')
    # get raw meta data for statistic
    meta <- SeuratS4@meta.data[,c("Sample","percent.mt", "percent.ribo",
                                  "nFeature_RNA","nCount_RNA")]
    # filter and normalize for detecting doublet.
    SeuratS4 <- subset(SeuratS4, subset = nFeature_RNA > param[['min.features']] & 
                         nFeature_RNA < param[['max.features']] & 
                         nCount_RNA > param[['min.UMIs']] &
                         nCount_RNA < param[['max.UMIs']] &
                         percent.mt < param[['max.mt']])
    SeuratS4 <- NormalizeData(SeuratS4,
                              normalization.method = "LogNormalize",
                              scale.factor = param[["scale.factor"]])
    # add noise gene score
    noiselist <- as.list(blacklist)
    noiselist <- lapply(noiselist, function(x) x[x != ''])
    allgene <- as.vector(unlist(noiselist))
    noiselist[['noise_gene']] <- allgene
    geneset <- c("Mitochondria", "Heat.shock.protein", "Ribosome", "Dissociation","allgenes")
    select_set <- 'noise_gene'
    regress_gene <- noiselist[select_set]
    SeuratS4 <- AddModuleScore(SeuratS4, features = regress_gene,
                               name = 'noise_gene')
    colnames(SeuratS4@meta.data)[colnames(SeuratS4@meta.data) == 'noise_gene1'] <- 'noise_gene'
    # hvg
    SeuratS4 <- FindVariableFeatures(object = SeuratS4,
                                     selection.method = param[['vargene.method']],
                                     nfeatures = param[['nFeatures']])
    if ('NULL' %in% param[['var.to.regress']]){
      vars.to.regress = NULL
    }else{
      vars.to.regress <- param[['var.to.regress']]
    }
    SeuratS4 <- ScaleData(SeuratS4, vars.to.regress = vars.to.regress)
    SeuratS4 <- RunPCA(SeuratS4)
    SeuratS4 <- RunUMAP(SeuratS4, dims =  param[['select_PCs']])
    # Seurat to CellDataSet object
    #message('Convert Seurat object to monocle2 object...')
    #cds <- as.CellDataSet(SeuratS4, assay = 'RNA', reduction = 'umap')
    message("Running scrublet...")
    scrublet_res <- scrublet_R(SeuratS4, python_home = param[['python_home']], 
                               return_results_only = T,
                               min_cells = param[['min.cells']],
                               expected_doublet_rate = param[['expected.doublet.rate']])
    message("Scrublet done.")
    SeuratS4$doublet_scores <- scrublet_res$doublet_scores
    p <- ggplot(data.frame(doublet_scores = SeuratS4$doublet_scores), 
                aes(x = doublet_scores)) + geom_histogram(bins = 200) + 
      scale_y_log10() + theme_bw()
    p <- p + geom_vline(xintercept = param[['max.doublet.rate']],
                        linetype="dashed", color = "black", size=1)
    ggsave(paste0(outdir, 'hist_doublet.pdf'),p, width = 6, height = 3.2)
    dbpos <- SeuratS4$doublet_scores > param[['max.doublet.rate']]
    sum(dbpos)
    #[1] 1529
    SeuratS4$doublet <- 'singlet'
    SeuratS4$doublet[dbpos] <- 'doublet'
    p1 <- FeaturePlot2(SeuratS4, features = 'doublet_scores')
    p2 <- DimPlot(SeuratS4, group.by = 'doublet')
    SeuratMajorVersion <- gsub('\\..*$','',packageVersion('Seurat'))
    device <- ifelse(SeuratMajorVersion == 4, yes = 'pdf', no = 'png')
    ggsave(paste0(outdir,'doublet_score.',device), p1 + p2, device = device,
           width = 7.7, height = 3.9)
    saveRDS(SeuratS4, file = paste0(outdir,'rawObject.rds'))
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
    write.table(statmat, file = paste0(outdir, 'sample_cell_number.txt'), 
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
    SeuratS4 <- subset(SeuratS4, subset = nFeature_RNA > param[['min.features']] & 
                         nFeature_RNA < param[['max.features']] & 
                         nCount_RNA > param[['min.UMIs']] &
                         nCount_RNA < param[['max.UMIs']] &
                         percent.mt < param[['max.mt']])
    SeuratS4 <- NormalizeData(SeuratS4,
                              normalization.method = "LogNormalize",
                              scale.factor = param[["scale.factor"]])
    # add noise gene score
    noiselist <- as.list(blacklist)
    noiselist <- lapply(noiselist, function(x) x[x != ''])
    allgene <- as.vector(unlist(noiselist))
    noiselist[['noise_gene']] <- allgene
    geneset <- c("Mitochondria", "Heat.shock.protein", "Ribosome", "Dissociation","allgenes")
    select_set <- 'noise_gene'
    regress_gene <- noiselist[select_set]
    SeuratS4 <- AddModuleScore(SeuratS4, features = regress_gene,
                               name = 'noise_gene')
    colnames(SeuratS4@meta.data)[colnames(SeuratS4@meta.data) == 'noise_gene1'] <- 'noise_gene'
    #hvg
    SeuratS4 <- FindVariableFeatures(object = SeuratS4,
                                     selection.method = param[['vargene.method']],
                                     nfeatures = param[['nFeatures']])
    if ('NULL' %in% param[['var.to.regress']]){
      vars.to.regress = NULL
    }else{
      vars.to.regress <- param[['var.to.regress']]
    }
    SeuratS4 <- ScaleData(SeuratS4, vars.to.regress = vars.to.regress)
    SeuratS4 <- RunPCA(SeuratS4)
    set.seed(123)
    SeuratS4 <- RunUMAP(SeuratS4, dims =  param[['select_PCs']])
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
      sce <- scDblFinder::scDblFinder(GetAssayData(SeuratS4, slot="counts"), 
                         clusters = Idents(SeuratS4),
                         samples = SeuratS4$Sample,
                         BPPARAM = bp,
                         dbr = param[['scdblFinder_dbr']],
                         dbr.sd = param[['scdblFinder_dbr.sd']]
      )
    }else{
      set.seed(123)
      sce <- scDblFinder::scDblFinder(GetAssayData(SeuratS4, slot="counts"), 
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
    ggsave(paste0(outdir,'doublet_score.',device), p1 + p2, device = device,
           width = 7.7, height = 3.9)
    saveRDS(SeuratS4, file = paste0(outdir,'rawObject.rds'))
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
    write.table(statmat, file = paste0(outdir, 'sample_cell_number.txt'), 
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
    write.table(statmat, file = paste0(outdir, 'sample_cell_number.txt'), 
                quote = F, sep = '\t')
    print(dim(SeuratS4))
    SeuratS4 <- subset(SeuratS4, subset = nFeature_RNA > param[['min.features']] & 
                         nFeature_RNA < param[['max.features']] & 
                         nCount_RNA > param[['min.UMIs']] &
                         nCount_RNA < param[['max.UMIs']] &
                         percent.mt < param[['max.mt']])
  }
  
  #QC plot after filtering
  message('Running QC after filtering.')
  num_g <- length(unique(SeuratS4@meta.data$orig.ident))
  features1 = c('percent.mt','nFeature_RNA',"nCount_RNA",
                "percent.ribo", "percent.HSP", 
                "S.Score","G2M.Score")
  #p <- RidgePlot(SeuratS4, features = features1, ncol = 1,group.by = 'Sample') + ylab('')+
  #  theme(legend.position = 'none')
 # ggsave(paste(outdir, "QC.RidgePlot_density_filtered.pdf", sep = ""),p,
 #        width = 6, height = (8 + num_g * 1.2))
  
  #################################################
  # 5. Normalization
  message('5. Normalization.')
  if ('NULL' %in% param[['var.to.regress']]){
    vars.to.regress = NULL
  }else{
    vars.to.regress <- param[['var.to.regress']]
  }
  # add noise gene score
  noiselist <- as.list(blacklist)
  noiselist <- lapply(noiselist, function(x) x[x != ''])
  allgene <- as.vector(unlist(noiselist))
  noiselist[['noise_gene']] <- allgene
  geneset <- c("Mitochondria", "Heat.shock.protein", "Ribosome", "Dissociation","noise_gene")
  select_set <- 'noise_gene'
  regress_gene <- noiselist[select_set]
  SeuratS4 <- AddModuleScore(SeuratS4, features = regress_gene,
                             name = 'noise_gene')
  colnames(SeuratS4@meta.data)[colnames(SeuratS4@meta.data) == 'noise_gene1'] <- 'noise_gene'
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
      SeuratS4@assays$RNA@var.features <- SeuratS4@assays$RNA@var.features[!
        SeuratS4@assays$RNA@var.features %in% unlist(regress_gene)
      ]
    }
    if (sum(param[['var.to.regress']] %in% c('S.Score','G2M.Score')) > 1) {
      message('Excluding cellcycle genes in high variable genes.')
      SeuratS4@assays$RNA@var.features <- SeuratS4@assays$RNA@var.features[!
                                                                             SeuratS4@assays$RNA@var.features %in% unlist(cc.genes)
      ]
    }
    # Identify the 10 most highly variable genes
    top10 <- head(VariableFeatures(SeuratS4), 10)
    # plot variable features with and without labels
    plot1 <- VariableFeaturePlot(SeuratS4)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    pdf(paste(outdir,"Filter_genes.pdf",sep=''),width = 8,height = 7)
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
      scale.genes <- SeuratS4@assays$RNA@var.features
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
      SeuratS4@assays$SCT@var.features <- SeuratS4@assays$SCT@var.features[!SeuratS4@assays$SCT@var.features %in% unlist(regress_gene)]
    }
    if (sum(param[['var.to.regress']] %in% c('S.Score','G2M.Score')) > 1) {
      message('Excluding cellcycle genes in high variable genes.')
      SeuratS4@assays$SCT@var.features <- SeuratS4@assays$SCT@var.features[!
                                                                             SeuratS4@assays$SCT@var.features %in% unlist(cc.genes)
      ]
    }
  }else{
    stop('Only support normalization methods: LogNormalize or SCT.')
  }
  

  
  ###########################################
  # 8. Perform linear dimensional reduction (PCA)
  message('8. Perform linear dimensional reduction (PCA)')
  SeuratS4 <- RunPCA(SeuratS4, npcs = param[['npcs']], 
                     features = VariableFeatures(object = SeuratS4))
  pdf(paste(outdir,"Select_pc.pdf",sep=''),width = 8,height = 7)
  ElbowPlot(SeuratS4, ndims = param[['npcs']])
  dev.off()
  pdf(paste(outdir,"pca.pdf",sep=''),width = 8,height = 7)
  DimPlot(SeuratS4, reduction = "pca",group.by = 'Sample')
  dev.off()
  
  
  ##############################################
  # 9. Batch effect correction
  message('9. Batch effect correction.')
  # SeuratV3-CCA and harmony are available
  if (param[['batch_correct_method']] == 'SeuratV3') {
    # a. SeuratV3-CCA
    # split the dataset into a list of two seurat objects (stim and CTRL)
    objlist <- SplitObject(SeuratS4, split.by = param[['batch_correct_colname']])
    
    # normalize and identify variable features for each dataset independently
    objlist <- lapply(X = objlist, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = param[['vargene.method']],
                                nfeatures = param[['nFeatures']])
    })
    # select features that are repeatedly variable across datasets for integration
    features <- SelectIntegrationFeatures(object.list = objlist)
    obj_anchors <- FindIntegrationAnchors(object.list = objlist, anchor.features = features)
    # this command creates an 'integrated' data assay
    SeuratS4 <- IntegrateData(anchorset = obj_anchors)
    #Perform an integrated analysis
    
    # specify that we will perform downstream analysis on the corrected data note that the
    # original unmodified data still resides in the 'RNA' assay
    DefaultAssay(SeuratS4) <- "integrated"
    SeuratS4 <- ScaleData(SeuratS4, verbose = FALSE)
    SeuratS4 <- RunPCA(SeuratS4, npcs = param[['npcs']], 
                       features = VariableFeatures(object = SeuratS4),
                       verbose = FALSE)
    
  }
  if (param[['batch_correct_method']] == 'harmony') {
    # b. harmony
    print(paste0('Harmony version: ', packageVersion('harmony')))
    set.seed(123)
    SeuratS4 <- SeuratS4 %>% 
      harmony::RunHarmony(param[['batch_correct_colname']],  project.dim = F)
    param[['reduction_name']] <- 'harmony'
    
  }
  ###############################################
  # 10. Run Non-linear dimensional reduction
  message('10. Run Non-linear dimensional reduction.')
  set.seed(123)
  SeuratS4 <- RunUMAP(SeuratS4, reduction = param[['reduction_name']], 
                      dims =  param[['select_PCs']])
  SeuratS4@meta.data$UMAP_1 <- SeuratS4@reductions$umap@cell.embeddings[,'UMAP_1']
  SeuratS4@meta.data$UMAP_2 <- SeuratS4@reductions$umap@cell.embeddings[,'UMAP_2']
  
  #############################################
  # 11. Cluster the cells
  message('11. Cluster the cells.')
  set.seed(123)
  SeuratS4 <- FindNeighbors(SeuratS4, dims = param[['select_PCs']],
                            reduction = param[['reduction_name']], 
                            k.param = param[['k.param']])
  SeuratS4 <- FindClusters(SeuratS4, resolution = param[['resolution']],
                           algorithm = param[['cluster_algorithm']])
  p1 <- DimPlot(SeuratS4)
  p2 <- DimPlot(SeuratS4, group.by = 'Sample')
  ggsave(paste0(outdir, 'UMAP.pdf'), p1+p2, width = 12, height = 4.8)
  ############################################
  # 12. Do DEG analysis of raw clusters 
  if (param[['do_DEG_wilcox']] == 'TRUE') {
    message('12. Do DEG analysis of raw clusters.')
    mmk <- FindAllMarkers(SeuratS4, only.pos = T)
    write.table(mmk, quote = F, sep = '\t', 
                file = paste0(outdir, 'FindAllMarkers.xls'))
  }
  ############################################
  return(SeuratS4)
}

#' Parameter template for PHASE1_run_Seurat_v3_QC_clustering.
#' 
#' @export
PHASE1_run_Seurat_v3_QC_clustering_param_template <- function() {
  ##config
  param <- list()
  param[['samplelist']] <- "NULL"
  param[['object']] <- "./result/101.cluster_main/final.rds"
  param[['subset_colname']] <- 'maintypes'
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
  param[['cal_cellcycle']] <- 'TRUE'
  param[['plot.raw']] <- 'TRUE'
  param[['normalize_method']] <- 'LogNormalize'
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
  param[['batch_correct_method']] <- 'NULL'
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
  param[['detect.doublet']] <- 'NULL'
  param[['scdblFinder_dbr']] <- 'NULL'
  param[['scdblFinder_dbr.sd']] <- 'NULL'
  param[['filter.doublet']] <- F
  param[['python_home']] <- '/home/jasper/.conda/envs/easylook/bin/python'
  param[['expected.doublet.rate']] <- 0.05
  param[['max.doublet.rate']] <- 0.2
  param[['is_multidata']] <- 'TRUE'
  param[['species']] <- 'Human'
  param[['do_DEG_wilcox']] <- 'FALSE'
  param[['remove_noise']] <- 'FALSE'
  
  return(param)
}

