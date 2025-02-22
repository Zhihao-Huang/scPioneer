################## Annotation
#' Annotate raw clusters by max expression, and return an object with annotation Idents
#'
#' @examples
#' markerdf <- data.frame(celltypes = c('T','T','NK','NK','Mono', 'Mono','DC','DC','B','B','Platelet'),
#' markers = c('CD3D','CD3E','NCAM1','NKG7','CD14','FCGR3A','CST3','CD1C','CD79A','MS4A1','PPBP'))
#' colnames(markerdf)
#' #[1] "celltypes" "markers"
#' pbmc_anno <- anno_top_gene(pbmc, markerdf, group.by = 'seurat_clusters')
#' DimPlot(pbmc_anno, label = T)
#' @export
anno_top_gene <- function(object, markerdf, group.by = 'seurat_cluster'
) {
  df <- GetAssayData(object)
  df <- df[unique(markerdf$markers),]
  df <- t(as.matrix(df))
  ident <- as.vector(object@meta.data[,group.by])
  dflist <- list()
  for (i in unique(ident)) {
    avedf <- apply(df[ident == i,], 2, mean)
    dflist[[paste0('cluster_',i)]] <- avedf
  }
  data <- do.call(rbind, dflist)
  celldflist <- list()
  for (cellt in unique(markerdf$celltypes)) {
    genes <- unique(markerdf$markers[markerdf$celltypes == cellt])
    celldflist[[cellt]] <- apply(data[,genes, drop = F],1,mean)
  }
  celldf <- do.call(cbind, celldflist)
  scaledata <- scale(celldf)
  annotation <- sapply(rownames(scaledata), function(x) colnames(scaledata)[which.max(scaledata[x,])])
  annodf <- data.frame(seurat_clusters = gsub('cluster_','',names(annotation)),
                       Annotation = as.vector(annotation),
                       stringsAsFactors = F)
  colnames(annodf)[1] <- c(group.by)
  meta <- object@meta.data
  meta$Annotation <- NULL
  meta$Cellname <- rownames(meta)
  meta <- left_join(meta, annodf)
  rownames(meta) <- meta$Cellname
  meta$Annotation <- factor(meta$Annotation, levels = unique(markerdf$celltypes))
  Idents(object) <- meta$Annotation
  return(object)
}

#' Annotate Celltypes by SingleR
#'
#' @examples
#' pbmc <- anno_SingleR(pbmc, species = 'Human') # based on individual cells
#' pbmc <- anno_SingleR(pbmc, species = 'Human', raw_cluster = 'seurat_clusters') # based on seurat clusters
#' @export
anno_SingleR <- function(object, species, assay = 'RNA', raw_cluster = NULL,
                         ensembl_version = 98,
                         label_colname = 'SingleR_label_cell',
                         label_raw_cluster_colname = 'SingleR_label_raw_cluster',
                         ref = NULL,
                         ref.label = c("label.main", "label.fine", "label.ont"),
                         ...) {
  ref.label <-  match.arg(NULL, choices = ref.label)
  if (is.null(ref)) {
    hpca <- celldex::HumanPrimaryCellAtlasData(ensembl=F)
    bpe <- celldex::BlueprintEncodeData(ensembl=F)
    hpca$label.main <- paste0("HPCA.", hpca$label.main)
    bpe$label.main <- paste0("BPE.", bpe$label.main)
    shared <- intersect(rownames(hpca), rownames(bpe))
    ref <- SEtools::mergeSEs( list(se1=hpca[shared,], se2=bpe[shared,]) )
  }
  counts <- GetAssayData(object, assay = assay, layer = 'counts')
  if (species == "Mouse") {
    orthologs <- Get_orthologs_mouse_human(version = ensembl_version,#101 -> 98, 20230103
                                           remove.duplicated = F,
                                           only.one.to.one = T,
                                           using.local.data = T)
    enslist <- orthologs[!duplicated(orthologs$MGI_symbol),]
    rownames(enslist) <- enslist$MGI_symbol
    features <- rownames(counts)
    enslist <- enslist[features, ]
    counts <- counts[!is.na(enslist$HGNC_symbol),]
    enslist <- enslist[!is.na(enslist$HGNC_symbol),]
    pos <- duplicated(enslist$HGNC_symbol)
    counts <- counts[!pos,]
    enslist <- enslist[!pos,]
    rownames(counts) <- enslist$HGNC_symbol
  }
  sce <- Seurat::as.SingleCellExperiment(CreateSeuratObject(counts = counts))
  set.seed(1000)
  clusters <- scran::quickCluster(sce)
  sce.pbmc <- scran::computeSumFactors(sce, cluster=clusters)
  sce <- scater::logNormCounts(sce)
  pred <- SingleR::SingleR(test=sce, ref=ref, labels=ref@colData[,ref.label], ...)
  table(pred$labels)
  object@meta.data[,label_colname] <- pred$labels
  Idents(object) <- object@meta.data[,label_colname]
  if (!is.null(raw_cluster)) {
    meta <- object@meta.data[,c(label_colname, raw_cluster)]
    colnames(meta) <- c('label_colname', 'raw_cluster')
    id <- meta %>% group_by(raw_cluster) %>%
      summarize(label = names(which.max(table(label_colname))))
    cluster.ids.new <- id$label
    names(cluster.ids.new) <- id$raw_cluster
    Idents(object) <- object@meta.data[,raw_cluster]
    object <- RenameIdents(object, cluster.ids.new)
    object@meta.data[,label_raw_cluster_colname] <- Idents(object)
  }
  return(object)
}

#' Annotate Celltypes by AUCell
#'
#' @examples
#' pbmc <- anno_AUCell(pbmc, species = 'Human', scsig.subset = 'Bone_Marrow') # based on individual cells
#' pbmc <- anno_AUCell(pbmc, species = 'Human', , scsig.subset = 'Bone_Marrow', raw_cluster = 'seurat_clusters') # based on seurat clusters
#' @export
anno_AUCell <- function(object, species, assay = 'RNA', raw_cluster = NULL,
                        ensembl_version = 98,
                        label_colname = 'AUCell_label_cell',
                        label_raw_cluster_colname = 'AUCell_label_raw_cluster',
                        scsig = NULL,
                        scsig.subset = c('_','Cord_Blood','Esophagus','Stomach',
                        'Small_Intestine','Large_Intestine','PFC','Embryonic',
                        'Midbrain','Bone_Marrow','Liver','Fetal_Kidney','Adult_Kidney','Pancreas'),
                        nCores = 4,
                        ...) {
  if (is.null(scsig)) {
    bfc <- BiocFileCache::BiocFileCache(ask=FALSE)
    scsig.path <- BiocFileCache::bfcrpath(bfc, file.path("http://software.broadinstitute.org",
                                          "gsea/msigdb/supplemental/scsig.all.v1.0.symbols.gmt"))
    scsig <- GSEABase::getGmt(scsig.path)
  }
  counts <- GetAssayData(object, assay = assay, layer = 'counts')
  if (species == 'Mouse') {
    orthologs <- Get_orthologs_mouse_human(version = ensembl_version,#101 -> 98, 20230103
                                           remove.duplicated = F,
                                           only.one.to.one = T,
                                           using.local.data = T)
    enslist <- orthologs[!duplicated(orthologs$MGI_symbol),]
    rownames(enslist) <- enslist$MGI_symbol
    features <- rownames(counts)
    enslist <- enslist[features, ]
    counts <- counts[!is.na(enslist$HGNC_symbol),]
    enslist <- enslist[!is.na(enslist$HGNC_symbol),]
    pos <- duplicated(enslist$HGNC_symbol)
    counts <- counts[!pos,]
    enslist <- enslist[!pos,]
    rownames(counts) <- enslist$HGNC_symbol
  }
  rankings <- AUCell::AUCell_buildRankings(counts, nCores = nCores, plotStats=FALSE)

    # Restricting to the subset of scsig gene sets:
  scsig <- scsig[grep(scsig.subset, names(scsig))]

  # Applying MsigDB to the Muraro dataset, because it's human:
  scsig.aucs <- AUCell::AUCell_calcAUC(scsig, rankings, nCores = nCores, ...)
  scsig.results <- t(SummarizedExperiment::assay(scsig.aucs))
  full.labels <- colnames(scsig.results)[max.col(scsig.results)]
  table(full.labels)
  object@meta.data[,label_colname] <- full.labels
  Idents(object) <- object@meta.data[,label_colname]
  if (!is.null(raw_cluster)) {
    meta <- object@meta.data[,c(label_colname, raw_cluster)]
    colnames(meta) <- c('label_colname', 'raw_cluster')
    id <- meta %>% group_by(raw_cluster) %>%
      summarize(label = names(which.max(table(label_colname))))
    cluster.ids.new <- id$label
    names(cluster.ids.new) <- id$raw_cluster
    Idents(object) <- object@meta.data[,raw_cluster]
    object <- RenameIdents(object, cluster.ids.new)
    object@meta.data[,label_raw_cluster_colname] <- Idents(object)
  }
  return(object)
}

#' Annotate celltypes by ChatGPT.
#'
#' @examples
#' Idents(pbmc) <- pbmc$seurat_clusters
#' markers <- FindAllMarkers(pbmc, logfc.threshold = 0.5, test.use = 'MAST', only.pos = T)
#' top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
#' annodf <- anno_openai(top10)
#' @export
anno_openai <- function(deg = NULL, genelist = NULL, tissuename = NULL,
                        base_url = "http://chatapi.littlewheat.com/v1",
                        api_key = 'sk-HgtySiUAhSLiZTlDRhNE7aEbERJOuSumUveDxYfAUy8YvDfM',
                        model = "gpt-3.5-turbo",
                        seed = 123) {
  if (!is.null(deg)) {
    if (is.numeric(deg$cluster)) deg$cluster <- paste0('raw_cluster__',deg$cluster)
    all.cluster <- unique(deg$cluster)
    genelist <- lapply(all.cluster, function(x) deg$gene[deg$cluster == x])
    names(genelist) <- all.cluster
  }else if (is.null(genelist)) {
    stop('Please provide a data.frame from FindAllmarkers or a genelist.')
  }
  if (is.null(names(genelist))) names(genelist) <- paste0('raw_cluster__',length(genelist))
  input <- sapply(names(genelist), function(x) paste0(x, ':', paste(genelist[[x]], collapse = ','),'. '))
  content = paste(paste0("Identify cell types of ",
                         tissuename, " cells using the following markers separately for each"),
                  " row. Only provide the cell type name.",
                  " Some can be a mixture of multiple cell types.",
                  input, collapse = "\n")

  client <- openai::OpenAI( base_url = base_url,
                            api_key = api_key)
  completion <- client$chat$completions$create(
    model = model,
    seed = seed,
    messages = list(list("role" = "user", "content" = content))
  )
  text <- completion$choices[1][[1]]$message$content
  cat(text)
  text <-  gsub("(\n)\\1{0,}", "__celltype__", text)
  anno <- strsplit(text, split = '__celltype__')[[1]]
  anno <- gsub(' $','', anno)
  raw_clusters <- sapply(anno, function(x) gsub('\\:.*$|-.*$','',x))
  raw_clusters = gsub('raw_cluster__','',raw_clusters)
  anno_clusters <- sapply(anno, function(x) gsub('^.*\\:|^.*-|^.*\\.','',x))
  anno_clusters <- gsub('^ ','', anno_clusters)
  anno_clusters <- gsub('\\,.*$', '', anno_clusters)
  df <- data.frame(Orig_Idents = raw_clusters,
                   Celltype_predicted = anno_clusters)
  return(df)
}

#' Annotate cell types using an AI wrapper embedded in the ellmer R package.
#'
#' @examples
#' Idents(pbmc) <- pbmc$seurat_clusters
#' markers <- FindAllMarkers(pbmc, logfc.threshold = 0.5, test.use = 'MAST', only.pos = T)
#' top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
#' annodf <- anno_ellmer(top10)
#' # depends on local model ollama. First download and install Ollama, then install model with ollama pull llama3.2. Test on Mac.
#' annodf <- anno_ellmer(top10, llm_function = 'ollama')
#' @export
anno_ellmer <- function(deg = NULL, genelist = NULL, tissuename = NULL,
                        base_url = "http://chatapi.littlewheat.com/v1",
                        api_key = 'sk-HgtySiUAhSLiZTlDRhNE7aEbERJOuSumUveDxYfAUy8YvDfM',
                        model = "gpt-3.5-turbo",
                        llm_function = c('openai','deepseek','ollama'),
                        ollama_model = 'llama3.2',
                        prompts = NULL,
                        prompts.add = NULL,
                        rm_str = c('\\*'),
                        sep = c('\\: ','- ','\\* ','\\. '),
                        as.order = F,
                        return.prompt = F,
                        return.answer = F,
                        seed = 1234) {
  llm_function <- match.arg(NULL, choices = llm_function)
  if (!is.null(deg)) {
    if (is.numeric(deg$cluster)) deg$cluster <- paste0('raw_cluster__',deg$cluster)
    all.cluster <- unique(deg$cluster)
    genelist <- lapply(all.cluster, function(x) deg$gene[deg$cluster == x])
    names(genelist) <- all.cluster
  }else if (!is.null(genelist)) {
    all.cluster <- names(genelist)
  }else{
    stop('Please provide a data.frame from FindAllmarkers or a genelist.')
  }
  if (is.null(names(genelist))) names(genelist) <- paste0('raw_cluster__',length(genelist))
  input <- sapply(names(genelist), function(x) paste0(x, ':', paste(genelist[[x]], collapse = ','),'. '))
  if (is.null(prompts)) prompts <- c(list(paste0("Identify cell types of ",
                           tissuename, " cells using the following markers separately for each row."),
                    " Only provide one cell type name.",
                    'Do not interpret.',
                    " Some can be a mixture of multiple cell types."))
  content = paste(prompts, prompts.add, input, collapse = "\n")
  if (return.prompt) return(content)
  if (llm_function == 'openai') {
    chat <- ellmer::chat_openai(
      system_prompt = NULL,
      turns = NULL,
      base_url = base_url,
      api_key = api_key,
      model = model,
      seed = seed,
      api_args = list(),
      echo = c("none", "text", "all")
      )
    text <- chat$chat(content)
  }
  if (llm_function == 'deepseek') {
    chat <- ellmer::chat_openai(
      system_prompt = NULL,
      turns = NULL,
      base_url = base_url,
      api_key = api_key,
      model = model,
      seed = seed,
      api_args = list(),
      echo = c("none", "text", "all")
    )
    text <- chat$chat(content)
  }
  if (llm_function == 'ollama') {
      chat <- ellmer::chat_ollama(
        system_prompt = NULL,
        turns = NULL,
        base_url = "http://localhost:11434",
        model = ollama_model,
        seed = seed,
        api_args = list(),
        echo = NULL
      )
      # content <- paste0("Identify cell types of ",
      #                   tissuename, " cells using the following markers separately for each row.",
      # " Only provide one cell type name. Do not interpret.",
      # " Some can be a mixture of multiple cell types.")
      # for (i in input) {
      #   text <- chat$chat(paste0(content))
      # }
      text <- chat$chat(content)
  }
  if (return.answer) return(text)
  text <-  gsub(rm_str, "", text)
  text <-  gsub("(\n)\\1{0,}", "__celltype__", text)
  anno <- strsplit(text, split = '__celltype__')[[1]]
  anno <- gsub(' $','', anno)
  if (as.order) {
    raw_clusters <- all.cluster
  }else{
    raw_clusters <- sapply(anno, function(x) gsub(paste0(paste(sep, collapse = '.*$|'),'.*$'),'',x))
  }
  #print(raw_clusters)
  raw_clusters = gsub('^.*raw_cluster__','',raw_clusters)
  anno_clusters <- sapply(anno, function(x) gsub(paste0('^.*', paste(sep, collapse = '|^.*')),'',x))
  anno_clusters <- gsub('^ ','', anno_clusters)
  anno_clusters <- gsub('\\,.*$', '', anno_clusters)
  #print(anno_clusters)

  df <- data.frame(Orig_Idents = raw_clusters,
                   Celltype_predicted = anno_clusters)
  print(head(df))
  return(df)
}

#' Annotate cell types using LLM model and DEGs
#'
#'
#' @export
anno_llm <- function(object, deg, raw_cluster, label_raw_cluster_colname, llm_function,
                     base_url, api_key, ...) {
  if (!is.null(raw_cluster)) Idents(object) <- object@meta.data[, raw_cluster]
  df <- anno_ellmer(deg = deg, llm_function = llm_function, base_url = base_url,
                    api_key = api_key, ...)
  new.id <- df$Celltype_predicted
  names(new.id) <- df$Orig_Idents
  Idents(object) <- object@meta.data[,raw_cluster]
  object <- RenameIdents(object, new.id)
  object@meta.data[,label_raw_cluster_colname] <- Idents(object)
  return(object)
}


#' Annotate Celltypes in one function
#'
#' @examples
#' ### SingleR
#' pbmc <- FindClusters(pbmc, resolution = 3)
#' obj <- annocell(pbmc, species = 'Human', method = 'SingleR')
#' obj <- annocell(pbmc, species = 'Human', method = 'SingleR', raw_cluster = 'seurat_clusters')
#' DimPlot(obj, label = T)
#'
#' ### AUCell
#' obj <- annocell(pbmc, species = 'Human', method = 'AUCell', scsig.subset = 'Bone_Marrow')
#' obj <- annocell(pbmc, species = 'Human', method = 'AUCell', scsig.subset = 'Bone_Marrow', raw_cluster = 'seurat_clusters')
#' DimPlot(obj, label = T)
#'
#' ### top_gene
#' markerdf <- data.frame(celltypes = c('T','T','NK','NK','Mono', 'Mono','DC','DC','B','B','Platelet'),
#' markers = c('CD3D','CD3E','NCAM1','NKG7','CD14','FCGR3A','CST3','CD1C','CD79A','MS4A1','PPBP'))
#' colnames(markerdf)
#' # "celltypes" "markers"
#' obj <- annocell(pbmc, species = 'Human', method = 'topgene', markerdf = markerdf, raw_cluster = 'seurat_clusters')
#' DimPlot(obj, label = T)
#'
#' ### angrycell
#' obj <- annocell(pbmc, species = 'Human', method = 'angrycell', raw_cluster = 'seurat_clusters')
#' DimPlot(obj, label = T)
#'
#' ### angrycell & OpenAI
#' Idents(pbmc) <- pbmc$seurat_clusters
#' markers <- FindAllMarkers(pbmc, logfc.threshold = 0.5, test.use = 'MAST', only.pos = T)
#' top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
#' obj <- annocell(pbmc, species = 'Human', method = 'angrycell', db = 'openai', DE = top10, raw_cluster = 'seurat_clusters',
#'  model = "gpt-3.5-turbo", seed = 123,
#'  base_url = "http://chatapi.littlewheat.com/v1",
#'  api_key = 'sk-HgtySiUAhSLiZTlDRhNE7aEbERJOuSumUveDxYfAUy8YvDfM')
#' DimPlot(obj, label = T)
#' @export
annocell <- function(object, species, assay = 'RNA', raw_cluster = NULL,
                     method = c('SingleR','AUCell','llm','topgene','angrycell'),
                     llm_function = c('openai','deepseek','ollama'),
                     ensembl_version = 98,
                     label_colname = NULL,
                     label_raw_cluster_colname = NULL,
                     ref = NULL,
                     ref.label = c("label.main", "label.fine", "label.ont"),
                     scsig = NULL,
                     scsig.subset = c('_','Cord_Blood','Esophagus','Stomach',
                                      'Small_Intestine','Large_Intestine','PFC','Embryonic',
                                      'Midbrain','Bone_Marrow','Liver','Fetal_Kidney','Adult_Kidney','Pancreas'),
                     markerdf = NULL,
                     DE = NULL,
                     min.pct = 0.3,
                     db = NULL,
                     n.cores = 4,
                     ...) {
  method <- match.arg(NULL, choices = method)
  llm_function <- match.arg(NULL, choices = llm_function)
  ref.label <-  match.arg(NULL, choices = ref.label)
  scsig.subset <-  match.arg(NULL, choices = scsig.subset)
  if (!is.null(raw_cluster)) {
    Idents(object) <- object@meta.data[, raw_cluster]
  }else{
    message('Idents(object) were used for raw clusters, as raw_cluster is NULL.')
    raw_cluster <- 'raw_cluster_temp'
    object$raw_cluster_temp <- Idents(object)
  }

  if (is.null(label_colname)) label_colname <- paste0(method, '_label_cell')
  if (is.null(label_raw_cluster_colname)) label_raw_cluster_colname <- paste0(method, '_label_raw_cluster')
  if (method == 'SingleR') object <- anno_SingleR(object, species, assay, raw_cluster,
                                               ensembl_version,
                                               label_colname, label_raw_cluster_colname,
                                               ref, ref.label,
                                               BPPARAM = BiocParallel::MulticoreParam(n.cores),
                                               ...)
  if (method == 'AUCell') object <- anno_AUCell(object, species, assay, raw_cluster,
                                             ensembl_version,
                                             label_colname, label_raw_cluster_colname,
                                             scsig, scsig.subset,
                                             nCores = n.cores, ...)
  if (method == 'llm') object <- anno_llm(object, DE, raw_cluster,label_raw_cluster_colname,
                                          llm_function, ...)
  if (method == 'topgene') {
    object <- anno_top_gene(object, markerdf, group.by = raw_cluster)
    object@meta.data[,label_raw_cluster_colname] <- Idents(object)
  }
  if (method == 'angrycell') {
    if (is.null(db)) db <- ifelse(species == 'Human', 'panglaodb_hs', 'panglaodb_mm')
    df <- angrycell(object, select.db = db, show.top = 1, min.pct = min.pct, core.num = n.cores,
                    active.assay = assay,...)
    new.id <- df$Celltype_predicted
    names(new.id) <- df$Orig_Idents
    Idents(object) <- object@meta.data[,raw_cluster]
    object <- RenameIdents(object, new.id)
    object@meta.data[,label_raw_cluster_colname] <- Idents(object)
  }
  return(object)
}



