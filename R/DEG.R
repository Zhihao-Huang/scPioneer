pos_pct <- function(x) {
  sum( x > 0) / length(x)
}
fold_change <- function(x, y, epsilon = 0.001){
  (mean(x) + epsilon) / (mean(y) + epsilon)
}
#' Find pair DEG.
#'@export
pairwisede <- function(seurat, ident.1=NULL, ident.2=NULL, group.by=NULL, slot.use = "data", assay = "RNA",
                       test.use="t"){
  if (is.null(group.by)) {
    groups <- Idents(seurat)
  } else {
    groups <- seurat@meta.data[[group.by]]
  }
  mask.1 <- groups == ident.1
  if (is.null(ident.2)) {
    mask.2 <- groups != ident.1
  } else {
    mask.2 <- groups == ident.2
  }
  expmat.1 <- t(as.matrix(slot(seurat@assays[[assay]], slot.use)[, mask.1]))
  expmat.2 <- t(as.matrix(slot(seurat@assays[[assay]], slot.use)[, mask.2]))
  emat.1 <- lapply(seq_len(ncol(expmat.1)), function(i) expmat.1[,i])
  emat.2 <- lapply(seq_len(ncol(expmat.1)), function(i) expmat.2[,i])
  if (test.use == "t"){
    stat <- mapply(t.test, emat.1, emat.2)
  } else if (test.use == "wilcox") {
    stat <- mapply(wilcox.test, emat.1, emat.2)
  }
  tibble(
    gene = rownames(slot(seurat@assays[[assay]], slot.use)),
    pct.1 = sapply(emat.1, pos_pct),
    pct.2 = sapply(emat.2, pos_pct),
    avg_logFC = log2(mapply(fold_change, emat.1, emat.2)),
    p_val = unlist(stat["p.value",]),
    p_val_adj = p.adjust(p_val, method = "fdr")
  )
}
#' Find different expression genes.
#'
#' @return A dataframe of markers list.
#' @examples
#' data('pbmc')
#' pbmc@@meta.data$Group <- Idents(pbmc)
#' markers <- allde(seurat=pbmc, group.by="Annotation", BPPARAM=MulticoreParam(12))
#' @export
allde <- function(seurat, group.by=NULL, BPPARAM=SerialParam(), ...){
  if (is.null(group.by)) {
    groups <- Idents(seurat)
  } else {
    groups <- seurat@meta.data[[group.by]]
  }
  gs <- unique(groups)
  df <- bplapply(gs, function(m){
    ddf <- pairwisede(seurat, ident.1 = m, group.by = group.by, ...)
    ddf$cluster <- m
    ddf
  }, BPPARAM = BPPARAM)
  names(df) <- gs
  df
}

signDE <- function(seurat, group.by=NULL, BPPARAM=SerialParam(), min.pct = 0.25, only.pos = T,
                   min.avg_logFC = 0.25,test.use = 'wilcox', ...){
  if (is.null(group.by)) {
    groups <- Idents(seurat)
  } else {
    groups <- seurat@meta.data[[group.by]]
  }
  gs <- unique(groups)
  df <- bplapply(gs, function(m){
    ddf <- pairwisede(seurat, ident.1 = m, group.by = group.by, test.use = test.use, ...)
    ddf$cluster <- m
    ddf <- ddf[ddf$p_val < 0.05,]
    ddf <- ddf[ddf$pct.1 > min.pct,]
    if(only.pos){
      ddf <- ddf[ddf$avg_logFC > min.avg_logFC,]
    }else{
      ddf <- ddf[abs(ddf$avg_logFC) > min.avg_logFC,]
    }
    ddf
  }, BPPARAM = BPPARAM)
  names(df) <- gs
  df
}



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

#' Pairwise DEG using limma-trend
#' # citation: https://ucdavis-bioinformatics-training.github.io/2019-single-cell-RNA-sequencing-Workshop-UCD_UCSF/scrnaseq_analysis/scRNA_Workshop-PART6.html
#' 
#' @param ... Parameter to GetAssayData.
#' @examples
#' DEG <- DEG_limma_trend(pbmc, 'Annotation', ident.1 = 'Memory CD4 T', ident.2 = 'Naive CD4 T')
#' # includes random effects
#' pbmc$batch <- c(rep('sample1',1000), rep('sample2',1000),rep('sample3',638))
#' DEG <- DEG_limma_trend(pbmc, 'Annotation', ident.1 = 'Memory CD4 T', ident.2 = 'Naive CD4 T', batch.colname = 'batch')
#' @export
DEG_limma_trend <- function(object, group.colname, ident.1, ident.2,  
                      batch.colname = NULL, min.cell = 1,
                      calcNormFactors = F, max.cpm.cut.off = 1, 
                      prior.count = 2,
                      duplicateCorrelation.name = NULL, 
                      trend = F,
                      sort.by = "P",...) {
  if (calcNormFactors) {
    counts <- as.matrix(GetAssayData(object, slot = 'counts', assay = 'RNA'))
    #Create DGEList object
    d0 <- DGEList(counts)
    #2. Preprocessing Calculate normalization factors
    d0 <- calcNormFactors(d0)
    #Filter low-expressed genes
    drop <- apply(cpm(d0), 1, max) < max.cpm.cut.off
    d <- d0[!drop,] 
    expr <- cpm(d, log=TRUE, prior.count = prior.count)
  }else{
    # Filter out genes that are less than min.cell for every cell in the object.
    expr <- as.matrix(GetAssayData(object,...))
    bad <- rowSums(expr) < min.cell
    expr <- expr[!bad,]
  }
  object$group.colname <- object@meta.data[,group.colname]
  object$group.colname <- as.vector(object$group.colname)
  object$group.colname <- check_replace_name(object$group.colname,
                                               '\\+| ',  '.', fixed = F)
  ident.1.adj <- check_replace_name(ident.1, '\\+| ',  '.', fixed = F)
  ident.2.adj <- check_replace_name(ident.2, '\\+| ',  '.', fixed = F)
  if (is.null(batch.colname)) {
    mm <- model.matrix(~0 + group.colname, data = object@meta.data)
  }else{
    object$batch.colname <- object@meta.data[,batch.colname]
    object$batch.colname <- as.vector(object$batch.colname)
    object$batch.colname <- check_replace_name(object$batch.colname,
                                                 '\\+| ',  '.', fixed = F)
    mm <- model.matrix(~0 + group.colname + batch.colname, data = object@meta.data)
  }
  if (!is.null(duplicateCorrelation.name)) {
    object$dupC.name <- object@meta.data[,duplicateCorrelation.name]
    object$dupC.name <- as.vector(object$dupC.name)
    object$dupC.name <- check_replace_name(object$dupC.name,
                                               '\\+| ',  '.', fixed = F)
    corfit <- duplicateCorrelation(expr, mm, block=object$dupC.name)
    fit <- lmFit(expr, mm, block=object$dupC.name, correlation=corfit$consensus)
  }else{
    fit <- lmFit(expr, mm)
  }
  x <- paste0('group.colname', ident.1.adj,' - ', 'group.colname', ident.2.adj)
  contr <- makeContrasts(contrasts = x, levels = colnames(coef(fit)))
  DEG <- contrasts.fit(fit, contrasts = contr)
  DEG <- eBayes(DEG, trend=trend)
  top.table <- topTable(DEG, sort.by = sort.by, n = Inf)
  top.table$gene <- rownames(top.table)
  top.table$cluster <- ident.1
  top.table$cluster[top.table$logFC < 0] <- ident.2
  return(top.table)
}
#' Pairwise DEG using limma-voom
#' 
#' @param ... Parameter to GetAssayData.
#' @examples
#' DEG <- DEG_limma_voom(pbmc, 'Annotation', ident.1 = 'Memory CD4 T', ident.2 = 'Naive CD4 T')
#' # includes random effects
#' pbmc$batch <- c(rep('sample1',1000), rep('sample2',1000),rep('sample3',638))
#' DEG <- DEG_limma_voom(pbmc, 'Annotation', ident.1 = 'Memory CD4 T', ident.2 = 'Naive CD4 T', batch.colname = 'batch')
#' @export
DEG_limma_voom <- function(object, group.colname, ident.1, ident.2,  
                      batch.colname = NULL, max.cpm.cut.off = 1,
                      duplicateCorrelation.name = NULL, 
                      ...) {
  counts <- as.matrix(GetAssayData(object, slot = 'counts', assay = 'RNA'))
  #Create DGEList object
  d0 <- DGEList(counts)
  #2. Preprocessing Calculate normalization factors
  d0 <- calcNormFactors(d0)
  #Filter low-expressed genes
  drop <- apply(cpm(d0), 1, max) < max.cpm.cut.off
  DEGobject <- d0[!drop,] 
  
  object$group.colname <- object@meta.data[,group.colname]
  object$group.colname <- as.vector(object$group.colname)
  object$group.colname <- check_replace_name(object$group.colname,
                                             '\\+| ',  '.', fixed = F)
  ident.1.adj <- check_replace_name(ident.1, '\\+| ',  '.', fixed = F)
  ident.2.adj <- check_replace_name(ident.2, '\\+| ',  '.', fixed = F)
  if (is.null(batch.colname)) {
    mm <- model.matrix(~0 + group.colname, data = object@meta.data)
  }else{
    object$batch.colname <- object@meta.data[,batch.colname]
    object$batch.colname <- as.vector(object$batch.colname)
    object$batch.colname <- check_replace_name(object$batch.colname,
                                               '\\+| ',  '.', fixed = F)
    mm <- model.matrix(~0 + group.colname + batch.colname, data = object@meta.data)
  }
  if (!is.null(duplicateCorrelation.name)) {
    object$dupC.name <- object@meta.data[,duplicateCorrelation.name]
    object$dupC.name <- as.vector(object$dupC.name)
    object$dupC.name <- check_replace_name(object$dupC.name,
                                           '\\+| ',  '.', fixed = F)
    corfit <- duplicateCorrelation(expr, mm, block=object$dupC.name)
    y <- voom(DEGobject, mm, plot = T)
    fit <- lmFit(y, mm, block=object$dupC.name, correlation=corfit$consensus)
  }else{
    y <- voom(DEGobject, mm, plot = T)
    fit <- lmFit(y, mm) 
  }
   
  x <- paste0('group.colname', ident.1.adj,' - ', 'group.colname', ident.2.adj)
  contr <- makeContrasts(contrasts = x, levels = colnames(coef(fit)))
  DEG <- contrasts.fit(fit, contrasts = contr)
  DEG <- eBayes(DEG)
  top.table <- topTable(DEG, sort.by = "P", n = Inf)
  top.table$gene <- rownames(top.table)
  top.table$cluster <- ident.1
  top.table$cluster[top.table$logFC < 0] <- ident.2
  return(top.table)
}
#' Running DEG of all clusters using limma. For each cluster, 
#' cells are compared with cells from all other clusters.
#' 
#' @param ... Parameter to GetAssayData.
#' @examples
#' # includes random effects
#' pbmc$batch <- c(rep('sample1',1000), rep('sample2',1000),rep('sample3',638))
#' DEG <- DEG_All_limma(pbmc, 'Annotation', batch.colname = 'batch', n.cores = 8)
#' @export
DEG_All_limma <- function(object, group.colname, only.pos = T, 
                          logFC.threhold = 0.25, adj.P.Val = 1,
                          batch.colname = NULL, n.cores = 4, ...) {
  expr <- GetAssayData(object,...)
  # Filter out genes that are 0 for every cell in this cluster
  bad <- Matrix::rowSums(expr) == 0
  expr <- expr[!bad,]
  object$group.colname <- object@meta.data[,group.colname]
  object$group.colname <- as.vector(object$group.colname)
  object$group.colname.adj <- check_replace_name(object$group.colname,
                                             '\\+|-| ',  '.', fixed = F)
  alltype_name <- unique(object$group.colname)
  DEGlist <- parallel::mclapply(alltype_name, function(x) {
    message(paste0('Processing ', x, ' vs All other clusters...'))
    subs4  <- object
    subs4$group.colname[subs4$group.colname != x] <- 'Others'
    subs4$group.colname.adj <- check_replace_name(subs4$group.colname,
                                                   '\\+|-| ',  '.', fixed = F)
    ident.1.adj <- check_replace_name(x, '\\+|-| ',  '.', fixed = F)
    if (is.null(batch.colname)) {
      mm <- model.matrix(~0 + group.colname.adj, data = subs4@meta.data)
    }else{
      subs4$batch.colname <- subs4@meta.data[,batch.colname]
      subs4$batch.colname <- as.vector(subs4$batch.colname)
      subs4$batch.colname <- check_replace_name(subs4$batch.colname,
                                                 '\\+|-| ',  '.', fixed = F)
      mm <- model.matrix(~0 + group.colname.adj + batch.colname, data = subs4@meta.data)
    }
    if (dim(expr)[2] > 20000) {
      expr <- biospiper::bigdataTOmat(expr)
    }
    fit <- lmFit(expr, mm)  
    pair <- paste0('group.colname.adj', ident.1.adj,' - ', 'group.colname.adjOthers')
    contr <- makeContrasts(contrasts = pair, levels = colnames(coef(fit)))
    DEG <- contrasts.fit(fit, contrasts = contr)
    DEG <- eBayes(DEG)
    top.table <- topTable(DEG, sort.by = "P", n = Inf)
    top.table$gene <- rownames(top.table)
    top.table$cluster <- x
    if (only.pos) top.table <- top.table[top.table$logFC > 0,]
    top.table <- top.table[abs(top.table$logFC) > logFC.threhold & 
                             top.table$adj.P.Val < adj.P.Val,]
    message(paste0(x, ' was done.'))
    return(top.table)
  },mc.cores = n.cores)
  DEGdf <- do.call(rbind, DEGlist)
  return(DEGdf)
}

#' Pairwise DEG using MAST to account for random effect.
#' 
#' Firstly provide log2(x + 1) transformed values (either raw counts or TPMs)
#' directly to MAST without doing normalization (just transformation). 
#' The "ngeneson" parameter in the model will account for "cell size" and act as a sort of normalization. 
#' Refer to https://github.com/kdzimm/PseudoreplicationPaper/issues/2
#' 
#' @examples 
#' pbmc$batch <- c(rep('sample1',1000), rep('sample2',1000),rep('sample3',638))
#' subs4 <- subset(pbmc, Annotation %in% c('CD14+ Mono','FCGR3A+ Mono'))
#' MASTdf <- DEG_MAST_RE(subs4, 'Annotation', ident.1 = 'CD14+ Mono', ident.2 = 'FCGR3A+ Mono', batch.colname = 'batch',freq_expressed = 0.2)
#' MASTdf <- MASTdf[MASTdf$FDR < 0.05 & MASTdf$avg_logFC > 0.15,]
#' @export
DEG_MAST_RE <- function(object, group.colname, ident.1, ident.2,  
                        batch.colname = NULL, freq_expressed = 0.1,
                        adjust.p.method = 'fdr',
                        Adaptive_thres = F) {
  object <- subset(object, cells = colnames(object)[object@meta.data[,group.colname] %in% c(ident.1, ident.2)])
  log2counts <- log2(as.matrix(GetAssayData(object, slot = 'counts', assay = 'RNA')) + 1)
  coldata <- object@meta.data[,c(batch.colname, group.colname)]
  fData <- data.frame(gene=rownames(log2counts))
  sca <- suppressMessages(MAST::FromMatrix(exprsArray=log2counts, cData=coldata, fData=fData))
  sca <- sca[MAST::freq(sca)>0,]
  if (Adaptive_thres) {
    message('Filter genes using adaptive thersholding.')
    thres <- MAST::thresholdSCRNACountMatrix(assay(sca), nbins = 20, min_per_bin = 30)
    assays(sca, withDimnames = FALSE) <- list(thresh=thres$counts_threshold, tpm=assay(sca))
  }
  expressed_genes <- MAST::freq(sca) > freq_expressed
  sca <- sca[expressed_genes,]
  message(paste0(ident.1, ' vs ', ident.2, ' remains ',dim(sca)[1],' genes.'))
  # running zlm
  cdr2 <- colSums(SummarizedExperiment::assay(sca)>0)
  SummarizedExperiment::colData(sca)$ngeneson <- scale(cdr2)
  SummarizedExperiment::colData(sca)$Group <-
    factor(SummarizedExperiment::colData(sca)[,group.colname])
  SummarizedExperiment::colData(sca)$Sample <-
    factor(SummarizedExperiment::colData(sca)[,batch.colname])
  emptydf <- data.frame(gene = NA, avg_logFC = NA, ci.hi = NA, ci.lo = NA, 
                        z = NA, P.Value = NA, FDR = NA, 
                        stringsAsFactors = F, check.names = F)
  zlmCond <- tryCatch(
    {
      suppressMessages(MAST::zlm(~ ngeneson + Group + (1 | Sample),
                                 sca, method='glmer',ebayes = F,
                                 strictConvergence = FALSE))
      # The return value of `readLines()` is the actual value 
      # that will be returned in case there is no condition 
      # (e.g. warning or error). 
      # You don't need to state the return value via `return()` as code 
      # in the "try" part is not wrapped inside a function (unlike that
      # for the condition handlers for warnings and error below)
    },
    error=function(cond) {
      message(paste("Error occurred when running zlm: ", ident.1, ' vs ', ident.2))
      message("Here's the original error message:")
      message(cond)
      message("Return an empty data frame instead of error message.")
      # Choose a return value in case of error
      return(emptydf)
    },
    finally={
      message(paste0("Processed zlm: ", ident.1, ' vs ', ident.2, ' (',dim(sca)[1],' genes).'))
    }
  )    
  if (class(zlmCond) == 'data.frame') {
    return(zlmCond)
  }
  contrast_name <- paste0('Group', sort(c(ident.1, ident.2))[2])
  summaryCond <- suppressMessages(MAST::summary(zlmCond, doLRT=contrast_name))
  summaryDt <- summaryCond$datatable
  fcHurdle <- merge(summaryDt[summaryDt$contrast==contrast_name
                              & summaryDt$component=='logFC', c(1,7,5,6,8)],
                    summaryDt[summaryDt$contrast==contrast_name
                              & summaryDt$component=='H', c(1,4)],
                    by = 'primerid')
  
  fcHurdle <- stats::na.omit(as.data.frame(fcHurdle))
  colnames(fcHurdle)[1] <- 'gene'
  colnames(fcHurdle)[2] <- 'avg_logFC'
  colnames(fcHurdle)[6] <- 'P.Value'
  
  #fcHurdle[, FDR := p.adjust(`P.Value`,'fdr')]
  if (nrow(fcHurdle) > 1) fcHurdle$FDR <- p.adjust(fcHurdle$P.Value, adjust.p.method)
  fcHurdle <- fcHurdle[order(fcHurdle$avg_logFC, decreasing = T), ]
  return(fcHurdle)
}

#' Pairwise DEG using glmmTMB to account for random effect.
#' 
#' Refer to https://github.com/kdzimm/PseudoreplicationPaper/blob/master/Type_1_Error/Type%201%20-%20Tweedie%20GLMM.Rmd
#' 
#' @examples 
#' pbmc$batch <- c(rep('sample1',1000), rep('sample2',1000),rep('sample3',638))
#' subs4 <- subset(pbmc, Annotation %in% c('CD14+ Mono','FCGR3A+ Mono'))
#' MASTdf <- DEG_MAST_RE(subs4, 'Annotation', ident.1 = 'CD14+ Mono', ident.2 = 'FCGR3A+ Mono', batch.colname = 'batch',freq_expressed = 0.2)
#' MASTdf <- MASTdf[MASTdf$FDR < 0.05 & MASTdf$avg_logFC > 0.15,]
#' @export
DEG_glmmTMB <- function(object, group.colname, ident.1, ident.2,  
                     batch.colname = NULL, min.ave.counts = 3,
                     adjust.p.method = 'fdr') {
  object <- subset(object, cells = colnames(object)[object@meta.data[,group.colname] %in% c(ident.1, ident.2)])
  genecounts <- as.matrix(GetAssayData(object, slot = 'counts', assay = 'RNA'))
  genecounts <- genecounts[which(apply(genecounts, 1, mean) > min.ave.counts), ] + 1
  genecounts <- log(sweep(genecounts,2,apply(genecounts,2,mean),'/'))
  coldata <- object@meta.data[,c(batch.colname, group.colname)]
  colnames(coldata) <- c('DonorID','Status')
  coldata$Status <- as.factor(coldata$Status)
  coldata$DonorID <- as.factor(coldata$DonorID)
  genecounts <- t(genecounts[,rownames(coldata)])
  allcells <- cbind(coldata,genecounts)
  
  fitmixed <- lapply( (ncol(coldata) + 1) : ncol(allcells),
                      function(x){glmmTMB::glmmTMB(allcells[,x] ~ Status + (1 | DonorID),
                                                   data=allcells,
                                                   family=glmmTMB::tweedie,
                                                   ziformula= ~0)})
  summaries <- lapply(fitmixed, summary)
  contrast_name <- paste0('Status', sort(c(ident.1, ident.2))[2])
  summaryCond <- suppressMessages(MAST::summary(fitmixed, doLRT='StatusneCRSwNP'))
  conds <- lapply(summaries, function(x){stats::coef(x)$cond[2,]})
  df <- as.data.frame(do.call(rbind, conds))
  df$gene <- colnames(genecounts)
  df <- stats::na.omit(as.data.frame(df))
  colnames(df)[2] <- 'Std.Error'
  colnames(df)[3] <- 'z'
  colnames(df)[4] <- 'P.Value'
  df$avg_log2FC <- df$Estimate
  df$avg_log2FC[df$avg_log2FC > 0] <- log2( df$avg_log2FC[df$avg_log2FC > 0] + 1)
  df$avg_log2FC[df$avg_log2FC < 0] <- -log2( -df$avg_log2FC[df$avg_log2FC < 0] + 1)
  if (nrow(df) > 1) df$FDR <- p.adjust(df$P.Value, adjust.p.method)
  df <- df[order(df$avg_log2FC, decreasing = T), ]
  return(df)
}


#' Pairwise DEG of multi-groups per cell type using Seurat::FindMarkers and return DEGlist.
#' 
#' @examples 
#' pbmc$batch <- c(rep('sample1',1000), rep('sample2',1000),rep('sample3',638))
#' subs4 <- subset(pbmc, Annotation %in% c('NK','DC','B'))
#' DEGlist <- DEG_group(subs4, groupname = 'Annotation',test.use = 'limma-trend')
#' DEGlist <- DEG_group(subs4, groupname = 'Annotation',test.use = 'limma-trend',batch.colname = 'batch')
#' @export
DEG_group <- function (object, groupname,test.use = c("wilcox","t","MAST","DESeq2",
                                                      "limma-trend","limma-voom",
                                                      "MAST_RE","glmmTMB"),
                       min.pct = 0.25, logfc.threshold = 0.2, # only for wilcox
                       batch.colname = NULL,  # only for limma
                       min.cell.num = 3,
                       g.cores = 1,
                       save.path = NULL,
                       annotation = NULL,
                       ...) {
  test.use <- match.arg(arg = NULL, choices = test.use)
  if (!is.null(levels(object@meta.data[, groupname]))) {
    allgroup <- levels(object@meta.data[, groupname])
  }else{
    allgroup <- unique(object@meta.data[, groupname])
  }
  num <- table(object@meta.data[, groupname])
  if (any(num < min.cell.num)) {
    removeg <- names(num)[num < min.cell.num]
    message(paste0("Warning: remove Group with few number (<",min.cell.num,"): ", 
                   paste(removeg, collapse = ",")))
    allgroup <- names(num)[num >= min.cell.num]
    
    if (length(allgroup) < 2) {
      message(paste0("Warning: only one group is available, pass."))
      if (test.use %in% c("wilcox","t","MAST","DESeq2")) {
        emptydf <- data.frame(p_val = NA, avg_logFC = NA, 
                              pct.1 = NA, pct.2 = NA, adj.P.Val = NA, gene = NA, 
                              cluster = NA, stringsAsFactors = F)
      }else if (test.use == 'MAST_RE'){
        emptydf <- data.frame(gene = NA, avg_logFC = NA, ci.hi = NA, ci.lo = NA, 
                              z = NA, P.Value = NA,adj.P.Val = NA,  cluster = NA,
                              stringsAsFactors = F, check.names = F)
      }else if (test.use == 'glmmTMB') {
        emptydf <- data.frame(Estimate = NA, Std.Error = NA, z = NA,
                              P.Value  = NA,  gene = NA, avg_log2FC = NA, 
                              adj.P.Val = NA,  cluster = NA, stringsAsFactors = F)
      }
      else{
        emptydf <- data.frame(logFC = NA, AveExpr = NA, 
                              t = NA, P.Value = NA, adj.P.Val = NA, B = NA, 
                              gene = NA, cluster = NA, stringsAsFactors = F)
      }
      
      return(list(emptydf))
      # next
    }
    
    cellID <- colnames(object)[object@meta.data[, groupname] %in% 
                                 allgroup]
    object <- subset(object, cells = cellID)
  }
  # running DEG
  if (test.use == 'limma-trend') {
    DEG_G4list <- lapply(allgroup[-length(allgroup)], function(x) {
      iter <- which(allgroup %in% x) + 1
      vs_g <- allgroup[iter:length(allgroup)]
      DEGlist <- lapply(vs_g, function(y) {
        DEG <- DEG_limma_trend(object, ident.1 = x, ident.2 = y,
                         group.colname = groupname, 
                         batch.colname = batch.colname,
                         ...)
        DEG$cluster <- paste0(x, "_vs_", y)
        DEG
      })
      names(DEGlist) <- paste0(x, "_", vs_g)
      DEGlist
    })
    DEG_G4 <- unlist(DEG_G4list, recursive = FALSE)
  }else if (test.use == 'limma-voom') {
    DEG_G4list <- lapply(allgroup[-length(allgroup)], function(x) {
      iter <- which(allgroup %in% x) + 1
      vs_g <- allgroup[iter:length(allgroup)]
      DEGlist <- lapply(vs_g, function(y) {
        DEG <- DEG_limma_voom(object, ident.1 = x, ident.2 = y,
                         group.colname = groupname,
                         batch.colname = batch.colname,
                         ...)
        DEG$cluster <- paste0(x, "_vs_", y)
        DEG
      })
      names(DEGlist) <- paste0(x, "_", vs_g)
      DEGlist
    })
    DEG_G4 <- unlist(DEG_G4list, recursive = FALSE)
  }else if (test.use == 'MAST_RE') {
    DEG_G4list <- parallel::mclapply(allgroup[-length(allgroup)], function(x) {
      iter <- which(allgroup %in% x) + 1
      vs_g <- allgroup[iter:length(allgroup)]
      DEGlist <- lapply(vs_g, function(y) {
        DEG <- DEG_MAST_RE(object, ident.1 = x, ident.2 = y,
                              group.colname = groupname,
                              batch.colname = batch.colname,
                              ...)
        comparison <- paste(sort(c(x, y)), collapse = "_vs_")
        DEG$cluster <- comparison
        colnames(DEG)[7] <- 'adj.P.Val'
        print(head(DEG))
        if (!is.null(save.path)) {
          if (is.null(annotation)) {
            annotation <- gsub(' ','_',as.character(unique(object$Annotation)))
          }
          openxlsx::write.xlsx(DEG, file = paste0(save.path, '/DEG_group_',annotation,'_',comparison,'.xlsx'))
        } 
        return(DEG)
      })
      names(DEGlist) <- as.vector( sapply(DEGlist, function(x) unique(x$cluster)) )
      return(DEGlist)
    }, mc.cores = g.cores)
    DEG_G4 <- unlist(DEG_G4list, recursive = FALSE)
  }else if (test.use == 'glmmTMB') {
    DEG_G4list <- parallel::mclapply(allgroup[-length(allgroup)], function(x) {
      iter <- which(allgroup %in% x) + 1
      vs_g <- allgroup[iter:length(allgroup)]
      DEGlist <- lapply(vs_g, function(y) {
        DEG <- DEG_glmmTMB(object, ident.1 = x, ident.2 = y,
                           group.colname = groupname,
                           batch.colname = batch.colname,
                           ...)
        comparison <- paste(sort(c(x, y)), collapse = "_vs_")
        DEG$cluster <- comparison
        colnames(DEG)[7] <- 'adj.P.Val'
        print(head(DEG))
        if (!is.null(save.path)) {
          if (is.null(annotation)) {
            annotation <- gsub(' ','_',as.character(unique(object$Annotation)))
          }
          openxlsx::write.xlsx(DEG, file = paste0(save.path, '/DEG_group_',annotation,'_',comparison,'.xlsx'))
        } 
        return(DEG)
      })
      names(DEGlist) <- as.vector( sapply(DEGlist, function(x) unique(x$cluster)) )
      return(DEGlist)
    }, mc.cores = g.cores)
    DEG_G4 <- unlist(DEG_G4list, recursive = FALSE)
  }
  else{
    # Seurat inserted methods
    DEG_G4list <- lapply(allgroup[-length(allgroup)], function(x) {
      iter <- which(allgroup %in% x) + 1
      # x is target group, vs_g are other groups
      vs_g <- allgroup[iter:length(allgroup)]
      DEGlist <- lapply(vs_g, function(y) {
        DEG <- FindMarkers(object, ident.1 = x, ident.2 = y, 
                           group.by = groupname, test.use = test.use, min.pct = min.pct, 
                           logfc.threshold = logfc.threshold, verbose = F)
        DEG$gene <- rownames(DEG)
        DEG$cluster <- paste0(x, "_vs_", y)
        colnames(DEG)[5] <- 'adj.P.Val'
        DEG
      })
      names(DEGlist) <- paste0(x, "_", vs_g)
      DEGlist
    })
    DEG_G4 <- unlist(DEG_G4list, recursive = FALSE)
  }
  
  return(DEG_G4)
}

#' DEG of pairwise groups per cell type.
#' 
#' @examples 
#' samples <- c(rep(c('sample1','sample2','sample3'),(ncol(pbmc) -1)/3),'sample1')
#' pbmc@@meta.data$Sample <- samples
#' DEGlist <- DEG_group_per_celltype(pbmc, celltype = 'Annotation',group = 'Sample',mc.cores = 1)
#' openxlsx::write.xlsx(DEGlist, file = './DEG_groups_per_celltype.xlsx')
#' 
#' @export
DEG_group_per_celltype <- function (object, celltype, group, test.use = "wilcox", 
                                    min.pct = 0.25, logfc.threshold = 0.2, 
                                    adj.P.Val = 1,
                                    batch.colname = NULL,
                                    mc.cores = 8,...) {
  allcelltype <- as.vector(sort(unique(object@meta.data[, celltype])))
  DEG_cell_G4 <- mclapply(1:length(allcelltype), function(i) {
    cell <- allcelltype[i]
    cellID <- colnames(object)[object@meta.data[, celltype] == 
                                 cell]
    subobj <- subset(object, cells = cellID)
    message(paste0("Cell number of ", cell, ": ", dim(subobj)[2]))
    DEGlist <- DEG_group(subobj, groupname = group, min.pct = min.pct, 
                         logfc.threshold = logfc.threshold, test.use = test.use,
                         batch.colname = batch.colname,...)
    if (is.null(DEGlist)) {
      message("Warning: non DEG was found.")
      DEGlist <- list(DEGlist)
    }
    DEGmat <- do.call(rbind, DEGlist)
    cat(paste0("\n", cell, " was done.\n\n"))
    return(DEGmat)
  }, mc.cores = mc.cores)
  DEG_cell_G4 <- lapply(DEG_cell_G4, function(x) {
    if (class(x) == 'data.frame') {x[x$adj.P.Val < adj.P.Val,]}
    else{return(x)}})
  names(DEG_cell_G4) <- allcelltype
  return(DEG_cell_G4)
}
