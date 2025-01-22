#' Combine multi vectors in different lengths, filled with NA.
#'
#' @param rawlist A list contains multi vectors in which elements have names.
#' @return A dataframe.
#' @examples
#' a <- c(1,3,5)
#' b <- c(1,2,3,4)
#' c <- c(1,2,4,5)
#' names(a) <- c('a','c','e')
#' names(b) <- c('a','b','c','d')
#' names(c) <- c('a','b','c','e')
#' rawlist <- list(a,b,c)
#' data <- cbind_fill(rawlist)
#' @export
cbind_fill <- function(rawlist){
  rown <- unique(unlist(sapply(rawlist, names)))
  # If elements do not have names, given index as names (default).
  if (any(is.null(rown))) {
    n <- names(rawlist)
    rawlist <- sapply(rawlist, function(x) {
      names(x) <- seq(1:length(x))
      x
    })
    names(rawlist) <- n
    rown <- unique(unlist(sapply(rawlist, names)))
  }
  matlist <- lapply(rawlist, as.matrix)
  clist <- lapply(matlist, function (x){
    missr <- rown[! rown %in% rownames(x)]
    mat <- matrix(0, length(missr), ncol(x))
    rownames(mat) <- missr
    matc <- rbind(x, mat)
    matc <- matc[rown,]})
  data.frame(clist)
}

#' Estimate  two side cut-off values of input-vector data.
#'
#' @param rawlist A list contains multi vectors in which elements have names.
#' @return A list.
#' @examples
#' data("pbmc")
#' Gene_counts <- pbmc@meta.data$nFeature_RNA
#' a <- density(Gene_counts,adjust = 1)
#' genecf <- findpeak_boo(Gene_counts)
#' gene_counts_cut_off = c(genecf$vp_cutoff,genecf$pv_cutoff)
#' print('gene count cut off: ')
#' print(gene_counts_cut_off)
#' #[1]  167 4381
#' hist(Gene_counts,prob=TRUE,breaks=200)
#' lines(density(Gene_counts,adjust = 0.5), col="blue", lwd=2) # add a density estimate with defaults
#' abline(v=gene_counts_cut_off[1],col='red')
#' abline(v=gene_counts_cut_off[2],col='red')
#' @export
findpeak_boo <- function(v,frac=0.1,adjust=1){
  a <- density(v,adjust = adjust)
  x <- a$y
  y <- a$x
  series <- v
  r <- rle(x)
  pks <- which(rep(x = diff(sign(diff(c(-Inf, r$values, -Inf)))) == -2,
                   times = r$lengths))
  boos <- which(rep(x = diff(sign(diff(c(-Inf, r$values, -Inf)))) == 2,
                    times = r$lengths))
  n_pks <- length(pks)
  n_boos <- length(boos)
  ##find largest valley-to-peak and peak-to-valley
  boos0 <- c(1,boos)
  v_p <- x[pks] - x[boos0]
  p_v <- x[pks[-length(pks)]] - x[boos0[-1]]
  if(length(p_v) == 0){
    print("Erorr in getting cut-off value, get 0 instead")
    return(list(
      "pks"=0,
      "boos"=0,
      "n_pks"=0,
      "n_boos"=0,
      "vp_cutoff"=0,
      "pv_cutoff"=0
    ))
    break
  }
  max_vp <-which(v_p == max(v_p))
  max_pv <- which(p_v == max(p_v))
  ##the ranges of the 2 largest gaps
  range_maxvp <- c(boos0[max_vp],pks[max_vp])
  range_maxpv <- c(pks[-length(pks)][max_pv],boos0[-1][max_pv])
  seq_vp <- x[range_maxvp[1]:range_maxvp[2]]
  seq_pv <- x[range_maxpv[1]:range_maxpv[2]]
  ##set fraction to cut one-side of delta gap-value distribution
  gap_vp <- seq_vp[-1] - seq_vp[-length(seq_vp)]
  delta_gap_vp <- abs(gap_vp[-1] - gap_vp[-length(gap_vp)])
  pos <- which(cumsum(delta_gap_vp) < sum(delta_gap_vp) * frac)
  vp_cutoff <- range_maxvp[1] + length(pos)
  gap_pv <- seq_pv[-length(seq_pv)] - seq_pv[-1]
  delta_gap_pv <- abs(gap_pv[-length(gap_pv)] - gap_pv[-1])
  pos <- which(cumsum(delta_gap_pv) > sum(delta_gap_pv) * (1-frac))
  pv_cutoff <- range_maxpv[2] - length(pos)
  pv <- round((y[pv_cutoff]-y[1]) * max(series)/(max(y)-min(y)),0)
  vp <- round((y[vp_cutoff]-y[1]) * max(series)/(max(y)-min(y)),0)
  return(list(
    "pks"=round(y[pks],0),
    "boos"=round(y[boos],0),
    "n_pks"=n_pks,
    "n_boos"=n_boos,
    "vp_cutoff"=vp,
    "pv_cutoff"=pv
  ))
}

findpeak_cut <- function(x, y, series) {
  r <- rle(x)
  pks <- which(rep(x = diff(sign(diff(c(-Inf, r$values, -Inf)))) == -2,
                   times = r$lengths))
  boos <- which(rep(x = diff(sign(diff(c(-Inf, r$values, -Inf)))) ==
                      2, times = r$lengths))
  n_pks <- length(pks)
  n_boos <- length(boos)
  ## find largest valley-to-peak and peak-to-valley
  boos0 <- c(1, boos)
  v_p <- x[pks] - x[boos0]
  p_v <- x[pks[-length(pks)]] - x[boos0[-1]]
  max_vp <- which(v_p == max(v_p))
  max_pv <- which(p_v == max(p_v))
  ## the ranges of the 2 largest gaps
  range_maxvp <- c(boos0[max_vp], pks[max_vp])
  range_maxpv <- c(pks[-length(pks)][max_pv], boos0[-1][max_pv])
  seq_vp <- x[range_maxvp[1]:range_maxvp[2]]
  seq_pv <- x[range_maxpv[1]:range_maxpv[2]]
  ## set 10% to cut one-side of delta gap-value distribution
  gap_vp <- seq_vp[-1] - seq_vp[-length(seq_vp)]
  delta_gap_vp <- abs(gap_vp[-1] - gap_vp[-length(gap_vp)])
  pos <- which(cumsum(delta_gap_vp) < sum(delta_gap_vp) * 0.1)
  vp_cutoff <- range_maxvp[1] + length(pos)
  gap_pv <- seq_pv[-length(seq_pv)] - seq_pv[-1]
  delta_gap_pv <- abs(gap_pv[-length(gap_pv)] - gap_pv[-1])
  pos <- which(cumsum(delta_gap_pv) > sum(delta_gap_pv) * 0.9)
  pv_cutoff <- range_maxpv[2] - length(pos)
  pv <- round((y[pv_cutoff] - y[1]) * max(series)/(max(y) - min(y)),
              0)
  vp <- round((y[vp_cutoff] - y[1]) * max(series)/(max(y) - min(y)),
              0)
  return(list(pks = round(y[pks], 0), boos = round(y[boos], 0), n_pks = n_pks,
              n_boos = n_boos, vp_cutoff = vp, pv_cutoff = pv))
}


#' Statistc of cell number before and after cell filtering.
#'
#' Statistc of cell number before and after cell filtering. mt_cut_off: 25, 15, 10, 5%
#' @return A dataframe of cell number in different filtering condition; hist plot of cutoff.
#' @param orig_metadata A metadata must contain information: nFeature_RNA,
#' percent_mt, nCount_RNA and Group.
#' @param input_gmin Munimum gene number for filtering cell.
#' @examples
#' samples <- c(rep(c('sample1','sample2','sample3'),(ncol(pbmc) -1)/3),'sample1')
#' pbmc@@meta.data$Sample <- samples
#' base_stat <- data_stat(pbmc@@meta.data, 100)
#' @export
data_stat <- function(orig_metadata, input_gmin,do.plot = FALSE) {
  fstat <- orig_metadata
  dlist <- list()
  ## genecount
  Gene_counts <- fstat$nFeature_RNA
  a <- density(Gene_counts, adjust = 1)
  genecf <- findpeak_cut(a$y, a$x, Gene_counts)
  prd_gene_counts_cut_off = c(genecf$vp_cutoff, genecf$pv_cutoff)
  print("gene count cut off: ")
  print(prd_gene_counts_cut_off)
  ## mt
  percent_mt <- fstat$percent.mt
  a <- density(percent_mt, adjust = 1)
  mtcf <- findpeak_cut(a$y, a$x, percent_mt)
  prd_mt_cut_off = mtcf$pv_cutoff
  dlist[["percent.mt_mt20"]] <- tapply(fstat$percent.mt, fstat$Sample,
                                       mean)
  ## gene_counts_cut_off=c(119,6500);mt_cut_off=25
  print("percnet mt cut off: ")
  print(prd_mt_cut_off)
  ## plot
  if (do.plot) {
    hist(Gene_counts, prob = TRUE, breaks = 200)
    lines(density(Gene_counts, adjust = 0.5), col = "blue", lwd = 2)  # add a density estimate with defaults
    abline(v = prd_gene_counts_cut_off[1], col = "red")
    abline(v = prd_gene_counts_cut_off[2], col = "red")
    hist(percent_mt, prob = TRUE, breaks = 200)
    lines(density(percent_mt, adjust = 0.5), col = "blue", lwd = 2)  # add a density estimate with defaults
    abline(v = prd_mt_cut_off, col = "red")
  }
  gminlist <- c(200, prd_gene_counts_cut_off[1], input_gmin)
  for (gmin in gminlist) {
    dlist <- list()
    dlist[["num_orig"]] <- table(fstat$Sample)
    dlist[["nCount_orig"]] <- tapply(fstat$nCount_RNA, fstat$Sample,
                                     mean)
    dlist[["nFeature_orig"]] <- tapply(fstat$nFeature_RNA, fstat$Sample,
                                       mean)
    dlist[["percent.mt_orig"]] <- tapply(fstat$percent.mt, fstat$Sample,
                                         mean)
    fstat1 <- fstat[fstat$nFeature_RNA > gmin, ]
    dlist[[paste0("num_gcount", gmin)]] <- table(fstat1$Sample)
    dlist[[paste0("nCount_gcount", gmin)]] <- tapply(fstat1$nCount_RNA,
                                                     fstat1$Sample, mean)
    dlist[[paste0("nFeature_gcount", gmin)]] <- tapply(fstat1$nFeature_RNA,
                                                       fstat1$Sample, mean)
    dlist[[paste0("percent.mt_gcount", gmin)]] <- tapply(fstat1$percent.mt,
                                                         fstat1$Sample, mean)
    # gmax <- gene_counts_cut_off[2] fstat <- fstat[fstat$nFeature_RNA
    # <gmax,] dlist[[paste0('num_gcount',gmax)]] <- table(fstat$Samples)
    # dlist[[paste0('nCount_gcount',gmax)]] <- tapply(fstat$nCount_RNA,
    # fstat$Samples,mean) dlist[[paste0('nFeature_gcount',gmax)]] <-
    # tapply(fstat$nFeature_RNA, fstat$Samples,mean)
    # dlist[[paste0('percent.mt_gcount',gmax)]] <- tapply(fstat$percent.mt,
    # fstat$Samples,mean)
  }
  mt_cut_off_list = c(25, 15, 10, 5)
  for (mt_cut_off in mt_cut_off_list) {
    fstat2 <- fstat1[fstat1$percent.mt < mt_cut_off, ]
    dlist[[paste0("num_mt", mt_cut_off)]] <- table(fstat2$Sample)
    dlist[[paste0("nCount_mt", mt_cut_off)]] <- tapply(fstat2$nCount_RNA,
                                                       fstat2$Sample, mean)
    dlist[[paste0("nFeature_mt", mt_cut_off)]] <- tapply(fstat2$nFeature_RNA,
                                                         fstat2$Sample, mean)
    dlist[[paste0("percent.mt_mt", mt_cut_off)]] <- tapply(fstat2$percent.mt,
                                                           fstat2$Sample, mean)
  }
  mt_cut_off = prd_mt_cut_off
  fstat1 <- fstat1[fstat1$percent.mt < mt_cut_off, ]
  dlist[[paste0("num_predicted_mt", mt_cut_off)]] <- table(fstat1$Sample)
  dlist[[paste0("nCount_predicted_mt", mt_cut_off)]] <- tapply(fstat1$nCount_RNA,
                                                               fstat1$Sample, mean)
  dlist[[paste0("nFeature_predicted_mt", mt_cut_off)]] <- tapply(fstat1$nFeature_RNA,
                                                                 fstat1$Sample, mean)
  dlist[[paste0("percent.mt_predicted_mt", mt_cut_off)]] <- tapply(fstat1$percent.mt,
                                                                   fstat1$Sample, mean)
  statdata <- do.call(cbind, dlist)
  return(statdata)
}

FastExpMean <- function(x) {
  x <- exp(x) - 1
  rowm <- apply(x, 1, mean)
  rowml <- log(rowm + 1)
  return(rowml)
}
FastLogVMR <- function(x) {
  x <- exp(x) - 1
  rowm <- apply(x, 1, mean)
  nnZero <- apply(x, 1, function(x) sum(x!=0))
  rowv <- apply(x, 1, function(x) (x - mean(x))^2)
  rowv <- (rowv + (ncol(x) - nnZero) * rowm^2) / (ncol(x) - 1)
  rowv <- log(rowv/rowm)
  return(rowv)
}
FastLogVMR <- function(x) {
  x <- exp(x) - 1
  rowm <- apply(x, 1, mean)
  rowv <- apply(x, 1, var)
  rowv <- log(rowv/rowm)
  return(rowv)
}

#' Find variable genes
#'
#' @export
FindVariableFeatures1 <- function(object,
                                  selection.method = c('vst','mean.var.plot','dispersion')[1],
                                  loess.span = 0.3,
                                  clip.max = 'auto',
                                  num.bin = 20,
                                  binning.method = c('equal_width', 'equal_frequency')[1],
                                  nfeatures = 2000,
                                  mean.cutoff = c(0.1, 8),
                                  dispersion.cutoff = c(1, Inf),
                                  assay = 'RNA'
) {
  if (selection.method == 'vst') {
    ###1. vst method.
    #mean
    data <- object@assays[[assay]]@counts
    hvf.info <- data.frame(mean = Matrix::rowMeans(x = data))
    #variance
    #hvf.info$variance <- rowVars(as.matrix(data))
    #as.matrix may cause memory problem when used for large data.
    hvf.info$variance <- apply(data, 1, var)
    #variance.expected
    hvf.info$variance.expected <- 0
    hvf.info$variance.standardized <- 0
    not.const <- hvf.info$variance > 0
    fit <- loess(
      formula = log10(x = variance) ~ log10(x = mean),
      data = hvf.info[not.const, ],
      span = loess.span
    )
    hvf.info$variance.expected[not.const] <- 10 ^ fit$fitted
    #variance.standardized
    sdnor <- function(x, mean.x, var.ex, clip.max){
      nor_x <- (x - mean.x)/var.ex
      #when value in var.ex equal to 0 ,return na. Then convert na to 0.
      nor_x[is.na(nor_x)] <- 0
      nor_x[nor_x > clip.max] <- clip.max
      sdx <- apply(nor_x, 1, var)
      return(sdx)
    }
    if (clip.max == 'auto') {
      clip.max <- sqrt(x = ncol(x = data))
    }
    hvf.info$variance.standardized <- sdnor(data,
                                            hvf.info$mean,
                                            sqrt(hvf.info$variance.expected),
                                            clip.max)
    #hvf.info$vst.variable <- order(hvf.info$variance.standardized) <2001
    topgenes <- hvf.info %>% top_n(nfeatures, variance.standardized)
    hvf.info$variable <- rownames(hvf.info) %in% rownames(topgenes)
    colnames(x = hvf.info) <- paste0('vst.', colnames(x = hvf.info))
  }else{
    ###2. mean.var.plot method
    data <- object@assays[[assay]]@data
    mean.function = FastExpMean
    dispersion.function = FastLogVMR
    if (!inherits(x = mean.function, what = 'function')) {
      stop("'mean.function' must be a function")
    }
    if (!inherits(x = dispersion.function, what = 'function')) {
      stop("'dispersion.function' must be a function")
    }
    ##Calculate the variance to mean ratio (VMR) in non-logspace (return answer in log-space)
    feature.mean <- mean.function(data)
    feature.dispersion <- dispersion.function(data)
    names(x = feature.mean) <- names(x = feature.dispersion) <- rownames(x = data)
    feature.dispersion[is.na(x = feature.dispersion)] <- 0
    feature.mean[is.na(x = feature.mean)] <- 0
    data.x.breaks <- switch(
      EXPR = binning.method,
      'equal_width' = num.bin,
      'equal_frequency' = c(
        -1,
        quantile(
          x = feature.mean[feature.mean > 0],
          probs = seq.int(from = 0, to = 1, length.out = num.bin)
        )
      ),
      stop("Unknown binning method: ", binning.method)
    )
    data.x.bin <- cut(x = feature.mean, breaks = data.x.breaks)
    names(x = data.x.bin) <- names(x = feature.mean)
    mean.y <- tapply(X = feature.dispersion, INDEX = data.x.bin, FUN = mean)
    sd.y <- tapply(X = feature.dispersion, INDEX = data.x.bin, FUN = sd)
    feature.dispersion.scaled <- (feature.dispersion - mean.y[as.numeric(x = data.x.bin)]) /
      sd.y[as.numeric(x = data.x.bin)]
    names(x = feature.dispersion.scaled) <- names(x = feature.mean)
    hvf.info <- data.frame(feature.mean, feature.dispersion, feature.dispersion.scaled)
    rownames(x = hvf.info) <- rownames(x = data)
    colnames(x = hvf.info) <- paste0('mvp.', c('mean', 'dispersion', 'dispersion.scaled'))
    #cut off
    means.use <- (hvf.info[, 1] > mean.cutoff[1]) & (hvf.info[, 1] < mean.cutoff[2])
    dispersions.use <- (hvf.info[, 3] > dispersion.cutoff[1]) & (hvf.info[, 3] < dispersion.cutoff[2])
    topfeatures <- rownames(x = hvf.info)[which(x = means.use & dispersions.use)]
    hvf.info$mvp.variable <-  means.use & dispersions.use
    if (selection.method == 'dispersion') {
      ##3. dispersion method
      hvf.info <- hvf.info[order(hvf.info$mvp.dispersion, decreasing = TRUE), , drop = FALSE]
      hvf.info$disp.variable <- rownames(x = hvf.info) %in% head(x = rownames(x = hvf.info), n = nfeatures)
    }
  }
  return(hvf.info)
}

#' Find variable genes
#' @examples
#' data <- object@@assays[[assay]]@@counts
#' hvg.info <- FindVariableFeatures2(pbmc,selection.method = 'vst')
#' data <- object@@assays[[assay]]@@data
#' hvg.info <- FindVariableFeatures2(pbmc,selection.method = 'mean.var.plot')
#' @export
FindVariableFeatures2 <- function(data,
                                  selection.method = c('vst','mean.var.plot','dispersion')[1],
                                  loess.span = 0.3,
                                  clip.max = 'auto',
                                  num.bin = 20,
                                  binning.method = c('equal_width', 'equal_frequency')[1],
                                  nfeatures = 2000,
                                  mean.cutoff = c(0.1, 8),
                                  dispersion.cutoff = c(1, Inf)
) {
  if (selection.method == 'vst') {
    ###1. vst method.
    ##mean
    hvf.info <- data.frame(mean = Matrix::rowMeans(x = data))
    ##variance
    #hvf.info$variance <- rowVars(as.matrix(data))
    #as.matrix may cause memory problem when used for large data.
    hvf.info$variance <- apply(data, 1, var)
    #variance.expected
    hvf.info$variance.expected <- 0
    hvf.info$variance.standardized <- 0
    not.const <- hvf.info$variance > 0
    fit <- loess(
      formula = log10(x = variance) ~ log10(x = mean),
      data = hvf.info[not.const, ],
      span = loess.span
    )
    hvf.info$variance.expected[not.const] <- 10 ^ fit$fitted
    #variance.standardized
    sdnor <- function(x, mean.x, var.ex, clip.max){
      nor_x <- (x - mean.x)/var.ex
      #when value in var.ex equal to 0 ,return na. Then convert na to 0.
      nor_x[is.na(nor_x)] <- 0
      nor_x[nor_x > clip.max] <- clip.max
      sdx <- apply(nor_x, 1, var)
      return(sdx)
    }
    if (clip.max == 'auto') {
      clip.max <- sqrt(x = ncol(x = data))
    }
    hvf.info$variance.standardized <- sdnor(data,
                                            hvf.info$mean,
                                            sqrt(hvf.info$variance.expected),
                                            clip.max)
    #hvf.info$vst.variable <- order(hvf.info$variance.standardized) <2001
    topgenes <- hvf.info %>% top_n(nfeatures, variance.standardized)
    hvf.info$variable <- rownames(hvf.info) %in% rownames(topgenes)
    colnames(x = hvf.info) <- paste0('vst.', colnames(x = hvf.info))
  }else{
    ###2. mean.var.plot method
    mean.function = FastExpMean
    dispersion.function = FastLogVMR
    if (!inherits(x = mean.function, what = 'function')) {
      stop("'mean.function' must be a function")
    }
    if (!inherits(x = dispersion.function, what = 'function')) {
      stop("'dispersion.function' must be a function")
    }
    ##Calculate the variance to mean ratio (VMR) in non-logspace (return answer in log-space)
    feature.mean <- mean.function(data)
    feature.dispersion <- dispersion.function(data)
    names(x = feature.mean) <- names(x = feature.dispersion) <- rownames(x = data)
    feature.dispersion[is.na(x = feature.dispersion)] <- 0
    feature.mean[is.na(x = feature.mean)] <- 0
    data.x.breaks <- switch(
      EXPR = binning.method,
      'equal_width' = num.bin,
      'equal_frequency' = c(
        -1,
        quantile(
          x = feature.mean[feature.mean > 0],
          probs = seq.int(from = 0, to = 1, length.out = num.bin)
        )
      ),
      stop("Unknown binning method: ", binning.method)
    )
    data.x.bin <- cut(x = feature.mean, breaks = data.x.breaks)
    names(x = data.x.bin) <- names(x = feature.mean)
    mean.y <- tapply(X = feature.dispersion, INDEX = data.x.bin, FUN = mean)
    sd.y <- tapply(X = feature.dispersion, INDEX = data.x.bin, FUN = sd)
    feature.dispersion.scaled <- (feature.dispersion - mean.y[as.numeric(x = data.x.bin)]) /
      sd.y[as.numeric(x = data.x.bin)]
    names(x = feature.dispersion.scaled) <- names(x = feature.mean)
    hvf.info <- data.frame(feature.mean, feature.dispersion, feature.dispersion.scaled)
    rownames(x = hvf.info) <- rownames(x = data)
    colnames(x = hvf.info) <- paste0('mvp.', c('mean', 'dispersion', 'dispersion.scaled'))
    #cut off
    means.use <- (hvf.info[, 1] > mean.cutoff[1]) & (hvf.info[, 1] < mean.cutoff[2])
    dispersions.use <- (hvf.info[, 3] > dispersion.cutoff[1]) & (hvf.info[, 3] < dispersion.cutoff[2])
    topfeatures <- rownames(x = hvf.info)[which(x = means.use & dispersions.use)]
    hvf.info$mvp.variable <-  means.use & dispersions.use
    if (selection.method == 'dispersion') {
      ##3. dispersion method
      hvf.info <- hvf.info[order(hvf.info$mvp.dispersion, decreasing = TRUE), , drop = FALSE]
      hvf.info$disp.variable <- rownames(x = hvf.info) %in% head(x = rownames(x = hvf.info), n = nfeatures)
    }
  }
  return(hvf.info)
}

#' Auto-select numbers of PCs and resolution
#'
#'Find cut off value by top 80% - 95% accumulative PCstdev with gaps > 0.001.
#'
#' @param object A seurat object.
#' @param plotfile The filepath of ElbowPlot file.
#' @param pc.perc.low The minimum value of accumulative PCstdev.
#' @param pc.perc.up The maximum value of accumulative PCstdev
#' @param gap.value the value of gap between two adjacent PCstdev.
#' @param ndims Total Number of PCs to compute and store (100 by default)
#' @param return.pc.r Return PCs and resolution instead of the range of PCs.
#' @return A vector which contains lower PCs and upper PCs.
#' @examples
#' get_numPCs_CUM(pbmc, plotfile=NULL, ndims = 100, return.pc.r = T)
#' get_numPCs_CUM(pbmc)
#' @export
get_numPCs_CUM <- function(object, do.plot=F,ndims = 100,pc.perc.low = 80,
                           pc.perc.up = 95, gap.value = 0.001,
                           return.pc.r = F){
  stdev <- object@reductions$pca@stdev
  stdev <- (stdev - min(stdev)) / (max(stdev) - min(stdev))
  gaps <- stdev[1:(length(stdev)-1)] - stdev[2:length(stdev)]
  gaps <- c(0,gaps)
  pc_perc <- cumsum(stdev)/sum(stdev) *100
  pc_cutoff_low <- sum(pc_perc < pc.perc.low & gaps > gap.value)
  pc_cutoff_up <- sum(pc_perc < pc.perc.up & gaps > gap.value)
  pc_c <- c(pc_cutoff_low,pc_cutoff_up)
  if (do.plot) {
    #pdf(plotfile,6,3.5)
    #plot
    p <- ElbowPlot(object,ndims = ndims)
    p <- p + geom_vline(xintercept = pc_c , linetype="dotted",
                        color = "red", size=1)
    print(p)
    #dev.off()
  }
  if (return.pc.r) {
    #select median PC as final PC
    PCs_num <-ceiling(median(seq(pc_c[1],pc_c[2]))[1])
    r <- round(cumsum(object@reductions$pca@stdev)[PCs_num]/60,1)
    return(c(PCs_num,r))
  }else{
    return(pc_c)
  }
}

#' Auto-select numbers of PCs and resolution in multi-sample data.
#' 
#' Find cut off PC values by top 80% - 95% cumsum of PCstdev with gaps > 0.001,
#' and optimize the result by Shannon Entopy. The meta.data of object must contain "Sample".
#' 
#' @param object A seurat object.
#' @param plotfile the filepath of ElbowPlot file.
#' @param ndims Total Number of PCs to compute and store (100 by default)
#' @param k.param Defines k for the k-nearest neighbor algorithm.
#' @param algorithm Parameter in FindClusters of Seurat.
#' @param select.PCs Input ranges of PCs for searching.
#' @param select_resolution Input ranges of resolution for searching. for example: c(0.5, 2.5).
#' @param reduction PC, tsne or umap.
#' @param normalize.entropy Normalize Shannon Entropy of cell type fraction
#' in each cluster.
#' @param plot.tree If TURE, Creates a plot of a clustering tree showing
#' the relationship between clusterings at different resolutions.
#' @param ... Deliver parameters to get_numPCs_CUM.
#' @param return.list Return Shannon Entopy list instead of PCs and resolution.
#' 
#' @return A vector which contains predicted PCs and resolution.
#' @examples
#' samples <- c(rep(c('sample1','sample2','sample3'),(ncol(pbmc) -1)/3),'sample1')
#' pbmc@@meta.data$Sample <- samples
#' SE_evaluation(pbmc)
#' SE_evaluation(pbmc,normalize.entropy = T)
#' SE_evaluation(pbmc,select.PCs = c(8,15))
#' SE_evaluation(pbmc,select.PCs = c(8,15),normalize.entropy = T)
#' @export
SE_evaluation <- function(object, k.param = 30, nn.eps = 0, ndims = 100,
                          algorithm = 1, select.PCs = NULL, normalize.entropy = F,
                          plot.tree = F, return.list = F, select_resolution = NULL,
                          reduction = 'PC', ...){
  if (is.null(select.PCs)) {
    pc_range <- get_numPCs_CUM(object,ndims = ndims,...)
    print(paste0('Searching the best top PC between ',pc_range[1],' and ',pc_range[2],'...'))
  }else{
    pc_range <- select.PCs
    print(paste0('Searching the best top PC between ',pc_range[1],' and ',pc_range[2],'...'))
  }
  shannonE_list <- list()
  pc_low <- pc_range[1]
  pc_up <- pc_range[2]
  for (PCs_num in pc_low : pc_up){
    object <- FindNeighbors(object, dims = 1:PCs_num, k.param = k.param, nn.eps = nn.eps,verbose = F)
    if (is.null(select_resolution)) {
      initial_r <- round(cumsum(object@reductions$pca@stdev)[PCs_num]/60,1)
      if (initial_r < 0.5) {initial_r <- 0.5}
      #set a range of resolution: initial_r - 0.5 ~ initial_r + 0.5
      r_range <- cumsum(rep(0.1,10)) + (initial_r - 0.5)
    }else{
      r_range <- seq(select_resolution[1], select_resolution[2], 0.1)
    }
    object <- FindClusters(object, resolution = r_range, algorithm = algorithm,verbose = F)
    meta <- object@meta.data
    meta$Annotation <- as.vector(meta$Annotation)
    if (reduction == 'PC') {
      meta$PC1 <- object@reductions$pca@cell.embeddings[,'PC_1']
      meta$PC2 <- object@reductions$pca@cell.embeddings[,'PC_2']
      x_value = 'PC1'
      y_value = 'PC2'
    }
    if (reduction == 'tsne') {
      meta$tSNE_1 <- object@reductions$tsne@cell.embeddings[,'tSNE_1']
      meta$tSNE_2 <- object@reductions$tsne@cell.embeddings[,'tSNE_2']
      x_value = 'tSNE_1'
      y_value = 'tSNE_2'
    }
    if (reduction == 'umap') {
      meta$UMAP_1 <- object@reductions$umap@cell.embeddings[,'UMAP_1']
      meta$UMAP_2 <- object@reductions$umap@cell.embeddings[,'UMAP_2']
      x_value = 'UMAP_1'
      y_value = 'UMAP_2'
    }
    ##plot rosolution tree
    if (plot.tree) {
      p <- clustree(meta, prefix = "RNA_snn_res.")
      p1 <- clustree_overlay(meta, prefix = "RNA_snn_res.", x_value = x_value, y_value = y_value)
      print(p)
      print(p1)
    }
    ##caculate Entropy for clusters of each resolution.
    for (r in r_range) {
      meta$seurat_clusters <- meta[,paste0('RNA_snn_res.',r)]
      frac <- meta %>% group_by(seurat_clusters,Sample) %>% summarise(count = n()) %>% mutate(count/sum(count))
      colnames(frac) <- c('seurat_clusters','Sample','Counts','Fraction')
      if (normalize.entropy) {
        shannonE_of_each_cluster <- tapply(frac$Fraction,list(frac$seurat_clusters),
                                           function(x){entropy::entropy(x)/length(x)})
      }else{
        #not normalize shannon entropy:
        shannonE_of_each_cluster <- tapply(frac$Fraction,list(frac$seurat_clusters),
                                           function(x){entropy::entropy(x)})
      }
      name <- paste0('PC_',PCs_num,'_',r)
      shannonE_list[[name]] <- shannonE_of_each_cluster
    }
  }
  #E <- round(mean(shannonE_list),2)
  if (return.list) {
    return(shannonE_list)
  }else{
    SEmat <- cbind_fill(shannonE_list)
    SEmean <- apply(SEmat,2,function(x){mean(x[x!=0])})
    #select the median index of all the same max values as result.
    PC_r <- names(SEmean)[median(which(SEmean %in% max(SEmean)))]
    PC_r <- unlist(strsplit(PC_r, '_'))
    return(PC_r[2:3])
  }
}

#' Auto-select numbers of PCs and resolution.
#'
#'find cut off value by P-value by JackStrawPlot.\
#'
#' @param object A seurat object.
#' @param plotfile the filepath of ElbowPlot file.
#' @param ndims Total Number of PCs to compute and store (100 by default)
#' @param p.val The cut off p.value from JackStraw.
#' @return A vector which contains predicted PCs and resolution.
#' @examples
#' get_numPCs_JA(pbmc,dim=10)
#' @export
get_numPCs_JA <- function(object,dim=100,p.val=0.001,num.replicate=100){
  if (nrow(object@reductions$pca@jackstraw$overall.p.values) == 0){
    object <- JackStraw(object, num.replicate = num.replicate, dim=dim)
    object <- ScoreJackStraw(object, dims = 1:dim)
  }
  pval <- object@reductions$pca@jackstraw$overall.p.values
  pc_select <- pval[,1][pval[,2] < p.val]
  PCs_num <- length(pc_select)
  initial_r = round(cumsum(object@reductions$pca@stdev)[PCs_num]/60,1)
  return(c(PCs_num,initial_r))
}

#' Etimate the Score of PCs accounted for genes.
#'
#' @param object A seurat object.
#' @param geneset Input gene set.Default is 'all',representing mt,ribosome,
#' and HSP.
#' @param species Species. Only Human and Mouse are permitted.
#' @param return.plot If ture, return ggplot object.
#' @param return.data If ture, return PCs x Score plotdata.
#' @return ggplot object or plotdata.
#' @examples
#' PCs_contribution(pbmc)
#' PCs_contribution(pbmc,geneset = c('CD1C','CD14'))
#' @export
PCs_contribution <- function(object,geneset = 'all',species = 'HomoSapiens',
                             return.plot = T,return.data = F){
  PCsmat <- object@reductions$pca@feature.loadings
  PCsmat <- abs(PCsmat)
  genelist <- list()
  if (species == 'HomoSapiens'){
    genelist[["mt"]] <- rownames(PCsmat)[grepl("^MT-",rownames(PCsmat))]
    genelist[["ribo"]] <- rownames(PCsmat)[grepl("^RP[SL]",rownames(PCsmat))]
    genelist[["HSP"]] <- rownames(PCsmat)[grepl("^HSP",rownames(PCsmat))]
  }else if(species == 'Mouse'){
    genelist[["mt"]] <- rownames(PCsmat)[grepl("mt-",rownames(PCsmat))]
    genelist[["ribo"]] <- rownames(PCsmat)[grepl("^Rp[sl]",rownames(PCsmat))]
    genelist[["HSP"]] <- rownames(PCsmat)[grepl("^Hsp",rownames(PCsmat))]
  }else{
    print('Only Human and Mouse are permitted!')
  }
  
  if ('all' %in% geneset) {
    if (length(genelist[["mt"]]) == 1){
      mat <- t(PCsmat[genelist[["mt"]],])
      mt <- apply(mat,2,sum)
    }else{
      mt <- apply(PCsmat[genelist[["mt"]],],2,sum)
    }
    if (length(genelist[["ribo"]]) == 1){
      mat <- t(PCsmat[genelist[["ribo"]],])
      ribo <- apply(mat,2,sum)
    }else{
      ribo <- apply(PCsmat[genelist[["ribo"]],],2,sum)
    }
    if (length(genelist[["HSP"]]) == 1){
      mat <- t(PCsmat[genelist[["HSP"]],])
      HSP <- apply(mat,2,sum)
    }else{
      HSP <- apply(PCsmat[genelist[["HSP"]],],2,sum)
    }
    data <- data.frame(mt,ribo,HSP,stringsAsFactors = F)
    data$PCs <- rownames(data)
    plotdata <- melt(data,id.vars = c('PCs'),measure.vars = colnames(data)[-ncol(data)])
    colnames(plotdata) <- c('PCs','Genes','Score')
  }else{
    pos <- geneset %in% rownames(PCsmat)
    if (sum(pos) == 0) {
      stop("ERROR: geneset must be 'all' or a vector of genesthat exist in High variable genes.")
    }else if (sum(!pos) != 0){
      print(paste0(sum(!pos),' Gene: ',paste0(geneset[!pos],collapse = ','),
                   ' not in High variable genes.'))
    }
    if (sum(pos) == 1){
      mat <- t(PCsmat[geneset[pos],])
      plotdata <- data.frame(apply(mat,2,sum),stringsAsFactors = F)
      colnames(plotdata) <- 'Score'
      plotdata$PCs <- rownames(plotdata)
    }else{
      plotdata <- data.frame(apply(PCsmat[geneset[pos],],2,sum),stringsAsFactors = F)
      colnames(plotdata) <- 'Score'
      plotdata$PCs <- rownames(plotdata)
    }
    
  }
  plotdata$PCs <- gsub('PC_','',plotdata$PCs)
  plotdata$PCs <- factor(plotdata$PCs,levels = sort(unique(as.numeric(plotdata$PCs))))
  if (return.plot & ('all' %in% geneset)) {
    p <- ggplot(plotdata,aes(PCs,Score,group = Genes)) + geom_line(aes(color=Genes))+
      theme_bw()+
      theme(panel.grid.major=element_line(colour=NA),
            panel.grid.minor=element_line(colour=NA))
  }else if (return.plot){
    p <- ggplot(plotdata,aes(PCs,Score,group = 1)) + geom_line()+
      theme_bw()+
      theme(panel.grid.major=element_line(colour=NA),
            panel.grid.minor=element_line(colour=NA))
  }
  if (return.data){
    return(plotdata)
  }else{
    return(p)
  }
}

#' Statistic of the filtered matrix.
#'
#' @param object A seurat object.
#' @param outdir The filedir to save statistic figures.
#' @return A Seurat object.
#' @examples
#' pbmc <- QC_filtered(pbmc,'./')
#' @export
QC_filtered <- function (SeuratS4, outdir, do_cellcycle = "TRUE") {
    if (packageVersion("scater") < "1.18.6") {
        warning("The package version of scater is too old, which supposed to be at least 1.18.6.")
    }
    mt_perc <- SeuratS4@meta.data$percent.mt
    ribo_perc <- SeuratS4@meta.data$percent.ribo
    HSP_perc <- SeuratS4@meta.data$percent.HSP
    num_g <- length(unique(SeuratS4@meta.data$Sample))
    pdf(paste0(outdir, "QC_violin_filtered.pdf"), (4 + num_g/5.5), 
        4)
    p1 <- VlnPlot(SeuratS4, features = c("nFeature_RNA", "nCount_RNA", 
        "percent.mt", "percent.ribo"), ncol = 2, pt.size = 0.05)
    print(p1)
    if (do_cellcycle == "TRUE") {
        p2 <- VlnPlot(SeuratS4, features = c("S.Score", "G2M.Score"))
        print(p2)
    }
    dev.off()
    pdf(paste(outdir, "QC.RidgePlot_density_filtered.pdf", sep = ""), 
        6, (4 + num_g * 1.2))
    p <- RidgePlot(SeuratS4, features = c("percent.mt", "percent.ribo", 
        "percent.HSP", "nFeature_RNA", "nCount_RNA"), ncol = 1)
    print(p)
    if (do_cellcycle == "TRUE") {
        p2 <- RidgePlot(SeuratS4, features = c("S.Score", "G2M.Score"))
        print(p2)
    }
    dev.off()
    pdf(paste(outdir, "QC.meta_density_filtered.pdf", sep = ""), 
        6, 3)
    p1 <- VlnPlot(SeuratS4, features = c("nFeature_RNA", "nCount_RNA", 
        "percent.mt", "percent.ribo"), ncol = 2, pt.size = 0.05)
    print(p1)
    if (do_cellcycle == "TRUE") {
        p2 <- VlnPlot(SeuratS4, features = c("S.Score", "G2M.Score"))
        print(p2)
    }
    dev.off()
    pdf(paste(outdir, "QC.RidgePlot_density_filtered.pdf", sep = ""), 
        6, (4 + num_g * 1.2))
    p <- RidgePlot(SeuratS4, features = c("percent.mt", "percent.ribo", 
        "percent.HSP", "nFeature_RNA", "nCount_RNA"), ncol = 1)
    print(p)
    if (do_cellcycle == "TRUE") {
        p2 <- RidgePlot(SeuratS4, features = c("S.Score", "G2M.Score"))
        print(p2)
    }
    dev.off()
    pdf(paste(outdir, "QC.meta_density_filtered.pdf", sep = ""), 
        6, 3)
    hist(mt_perc, prob = TRUE, breaks = 200)
    hist(ribo_perc, prob = TRUE, breaks = 200)
    hist(HSP_perc, prob = TRUE, breaks = 200)
    dev.off()
    pdf(paste(outdir, "QC.meta_scater_filtered.pdf", sep = ""), 
        8, 5)
    sce <- as.SingleCellExperiment(SeuratS4)
    sce <- scater::logNormCounts(sce)
    sce$Samples <- sce$ident
    if (do_cellcycle == "TRUE") {
        vars <- scater::getVarianceExplained(sce, variables = c("Samples", 
            "percent.mt", "percent.HSP", "nFeature_RNA", "nCount_RNA", 
            "S.Score", "G2M.Score"))
    }
    else {
        vars <- scater::getVarianceExplained(sce, variables = c("Samples", 
            "percent.mt", "percent.HSP", "nFeature_RNA", "nCount_RNA"))
    }
    print(scater::plotExplanatoryVariables(vars))
    dev.off()
    return(SeuratS4)
}

#' Statistic of the orignal matrix.
#'
#' @param object A seurat object.
#' @param outdir The filedir to save statistic figures.
#' @param species Species of the sample. Only Human and Mouse are permitted.
#' @param multi.samples If TRUE, data contains samples from different individual.
#' @return A Seurat object.
#' @examples
#' pbmc <- QC_raw(pbmc,'./')
#' @export
QC_raw <- function (object, outdir, species = "Human", do_cellcycle = "TRUE", 
                    plot.raw = 'TRUE', 
                    max.mt = 10, min.features = 200, max.features = Inf,
                    multi.samples = "TRUE") 
{
  if (species == "Human") {
    object[["percent.mt"]] <- PercentageFeatureSet(object, 
                                                   pattern = "^MT-")
    object[["percent.ribo"]] <- PercentageFeatureSet(object, 
                                                     pattern = "^RP[SL]")
    object[["percent.HSP"]] <- PercentageFeatureSet(object, 
                                                    pattern = "^HSP")
  }
  else if (species == "Mouse") {
    object[["percent.mt"]] <- PercentageFeatureSet(object, 
                                                   pattern = "mt-")
    object[["percent.ribo"]] <- PercentageFeatureSet(object, 
                                                     pattern = "^Rp[sl]")
    object[["percent.HSP"]] <- PercentageFeatureSet(object, 
                                                    pattern = "^Hsp")
  }
  else {
    stop("Only Human and Mouse are permitted!")
  }
  mtNA <- is.na(object[["percent.mt"]])
  rbNA <- is.na(object[["percent.ribo"]])
  hspNA <- is.na(object[["percent.HSP"]])
  if (any(mtNA) | any(rbNA) | any(hspNA)) {
    message("NA occurs percentage of MT/RIBO/HSP when nUMIs is 0, which was tranfer as 0 later.")
    object[["percent.mt"]][mtNA] <- 0
    object[["percent.ribo"]][rbNA] <- 0
    object[["percent.HSP"]][hspNA] <- 0
  }
  if (do_cellcycle == "TRUE") {
    object <- NormalizeData(object)
    if (species == "Human") {
      object <- CellCycleScoring(object = object, g2m.features = cc.genes$g2m.genes, 
                                 s.features = cc.genes$s.genes)
    }else{
      orthologs <- Get_orthologs_mouse_human(version = 98,#101 -> 98, 20230103
                                                        remove.duplicated = F,
                                                        only.one.to.one = T,
                                                        using.local.data = T)
      rownames(orthologs) <- orthologs$HGNC_symbol
      g2m.genes.mouse <- as.vector(unique(orthologs[cc.genes$g2m.genes,]$MGI_symbol))
      s.genes.mouse <- as.vector(unique(orthologs[cc.genes$s.genes,]$MGI_symbol))
      if (any(is.na(g2m.genes.mouse))) {
        message(paste0('Human G2M genes does not exist in mouse: ',
                       paste(cc.genes$g2m.genes[is.na(g2m.genes.mouse)],collapse = ',')))
        g2m.genes.mouse <- g2m.genes.mouse[!is.na(g2m.genes.mouse)]
      }
      if (any(is.na(s.genes.mouse))) {
        message(paste0('Human S genes does not exist in mouse: ',
                       paste(cc.genes$s.genes[is.na(s.genes.mouse)],collapse = ',')))
        s.genes.mouse <- s.genes.mouse[!is.na(s.genes.mouse)]
      }
      object <- CellCycleScoring(object = object, g2m.features = g2m.genes.mouse, 
                                 s.features = s.genes.mouse)
    }
  }
  rawlist <- list()
  if (plot.raw == 'TRUE') {
    message("================step1 statistic parameters================")
    parameters <- c("percent.mt", "percent.ribo", "percent.HSP", 
                    "nFeature_RNA", "nCount_RNA", "S.Score", "G2M.Score")
    parameters <- parameters[parameters %in% colnames(object@meta.data)]
    if (multi.samples == "TRUE") {
      Idents(object) <- object@meta.data$orig.ident
      features1 = parameters[-1]
      plist <- lapply(features1, function(x) {
        RidgePlot(object, features = x, ncol = 1,group.by = 'Sample') + ylab('')+
          theme(legend.position = 'none')
      })
      names(plist) <- features1
      plist[['percent.mt']] <- RidgePlot(object, features = 'percent.mt',
                                         ncol = 1, group.by = 'Sample') +
        geom_vline(xintercept = max.mt,
                   linetype="dotted", color = "black", size=1)+ ylab('')+
        theme(legend.position = 'none')
      plist[['nFeature_RNA']] <- RidgePlot(object, features = 'nFeature_RNA', 
                                           ncol = 1,group.by = 'Sample')+
        geom_vline(xintercept = min.features,
                   linetype="dotted", color = "black", size=1)+
        geom_vline(xintercept = max.features,
                   linetype="dotted", color = "black", size=1)+ ylab('')+
        theme(legend.position = 'none')
      features1 = parameters
      plist <- plist[features1]
      num_g <- length(unique(object@meta.data$Sample))
      ggsave(paste0(outdir, "QC.RidgePlot_density.pdf"), patchwork::wrap_plots(plist, ncol = 1),
          width = 6, height = (8 + num_g * 1.2))
      rawlist <- c(rawlist, plist)
    }
    message("================step2 plot density=================")
    pdf(paste(outdir, "QC.meta_density.pdf", sep = ""), 6, 3)
    hist(object@meta.data$percent.mt, prob = TRUE, breaks = 200)
    hist(object@meta.data$percent.ribo, prob = TRUE, breaks = 200)
    hist(object@meta.data$percent.HSP, prob = TRUE, breaks = 200)
    dev.off()
    df <- object@meta.data[,c('nFeature_RNA','percent.mt')]
    df$percent.NOTmt <- 100 - df$percent.mt
    p <- ggplot(df, aes(percent.NOTmt, log2(nFeature_RNA+1))) +
      geom_point(size = 0.1) +
      geom_density2d(aes(colour=..level..)) + 
      scale_colour_gradient(low="green",high="red") +
      geom_hline(yintercept = log2(min.features+1),
                 linetype="dashed", color = "black", size=1)+
      geom_hline(yintercept = log2(max.features+1),
                 linetype="dashed", color = "black", size=1)+
      geom_vline(xintercept = 100-max.mt,
                 linetype="dashed", color = "black", size=1)+
      theme_bw() + theme(legend.position = 'none')
    ggsave(paste0(outdir, "QC.scatter.png"),p,device = 'png', 
           width = 5.2, height = 4.6)
    df$log2nFeature_RNA <-  log2(df$nFeature_RNA+1)
    p <- ggPoint2(
      x = df$percent.NOTmt,
      y = df$log2nFeature_RNA,
      colorDensity = TRUE,
      continuousSet = "sambaNight",
      size = 0.1,
      xlabel = "percentage of Non-MT counts",
      ylabel = "log2(Number of genes + 1)"
    ) + geom_hline(yintercept = log2(min.features+1), lty = "dashed") +
      geom_vline(xintercept = 100-max.mt, lty = "dashed")
    ggsave(paste0(outdir, "QC.scatter2.png"),p,device = 'png', 
           width = 5.2, height = 4.6)
    rawlist[['QC.scatter2']] <- p
    message("================step3 plot violin================")
    num_g <- length(unique(object@meta.data$Sample)) 
    plot1 <- FeatureScatter(object, feature1 = "nCount_RNA", 
                            feature2 = "percent.mt")
    plot2 <- FeatureScatter(object, feature1 = "nCount_RNA", 
                            feature2 = "nFeature_RNA")
    p <- CombinePlots(plots = list(plot1, plot2))
    ggsave(paste0(outdir, "QC.violin.pdf"), p, width = (12 + num_g/5.5), 
          height = 4)
    rawlist[['Feature_mt_Scatter']] <- plot1
    rawlist[['Feature_Count_Scatter']] <- plot2
  }
  orig.metadata <- object@meta.data
  write.table(file = paste(outdir, "orig.metadata.txt", sep = ""), 
              orig.metadata, quote = F, sep = "\t")
  return(list(object = object, plotlist = rawlist))
}

#' QC. Add percent.ribo or percent.HSP to object；return object or Ridgeplot.
#' 
#' @examples 
#' Idents(pbmc) <- pbmc$orig.ident
#' QC_plot(pbmc, species = 'Human', add.line = T,minF = 200, maxF = 10000, percmt = 10)
#' pbmc <- subset(pbmc, subset = nFeature_RNA > 800)
#' QC_plot(pbmc, species = 'Human', add.line = F)
#' Idents(pbmc) <- pbmc$Annotation
#' QC_plot(pbmc,species = 'Human', add.line = F)
#' ##return object with percent.ribo or percent.HSP
#' pbmc <- QC_plot(pbmc, species = 'Human', return.obj = T)
#' @export

QC_plot <- function (SeuratS4, species = "Human", ident = NULL, add.line = F, 
    minF = 200, maxF = 10000, percmt = 10, color.use = NULL, 
    combine = T, return.obj = F) 
{
    message(species)
    if (species == "Human") {
        SeuratS4[["percent.mt"]] <- PercentageFeatureSet(SeuratS4, 
            pattern = "^MT-")
        SeuratS4[["percent.ribo"]] <- PercentageFeatureSet(SeuratS4, 
            pattern = "^RP[SL]")
        SeuratS4[["percent.HSP"]] <- PercentageFeatureSet(SeuratS4, 
            pattern = "^HSP")
        angry_species <- "HomoSapiens"
    }
    else if (species == "Mouse") {
        SeuratS4[["percent.mt"]] <- PercentageFeatureSet(SeuratS4, 
            pattern = "mt-")
        SeuratS4[["percent.ribo"]] <- PercentageFeatureSet(SeuratS4, 
            pattern = "^Rp[sl]")
        SeuratS4[["percent.HSP"]] <- PercentageFeatureSet(SeuratS4, 
            pattern = "^Hsp")
        angry_species <- "Mouse"
    }
    else {
        message("Only Human and Mouse are provided!")
    }
    if (return.obj) {
        return(SeuratS4)
    }
    if (is.null(ident)) {
        ident <- as.vector(unique((Idents(SeuratS4))))
    }
    if (is.null(color.use)) {
        color.use <- scPalette1(length(ident))
    }
    gglist <- RidgePlot(SeuratS4, features = c("percent.mt", 
        "percent.ribo", "nFeature_RNA", "nCount_RNA"), idents = ident, 
        ncol = 2, cols = color.use, combine = F)
    if (add.line) {
        gglist[[1]] <- gglist[[1]] + ylab("") + geom_vline(aes(xintercept = percmt), 
            linetype = "dashed")
        gglist[[2]] <- gglist[[2]] + ylab("")
        gglist[[3]] <- gglist[[3]] + ylab("") + geom_vline(aes(xintercept = minF), 
            linetype = "dashed")
        gglist[[3]] <- gglist[[3]] + ylab("") + geom_vline(aes(xintercept = maxF), 
            linetype = "dashed")
        gglist[[4]] <- gglist[[4]] + ylab("")
    }
    if (combine) {
        gglist[[1]] + gglist[[2]] + gglist[[3]] + gglist[[4]] + 
            patchwork::plot_layout(ncol = 2)
    }
    else {
        gglist
    }
}

#' Comparing number of cell before and after filtering.
#' 
#' @examples 
#' meta <- scPioneer::meta200[,c("Sample","percent.mt", "percent.ribo","nFeature_RNA","nCount_RNA")]
#' statmat <- stat_num(meta)
#' statmat <- stat_num(meta, min.param.list = list(nFeature_RNA = 200,nCount_RNA = 500),
#' max.param.list = list(percent.mt = 10, percent.ribo = 70))
#' 
#' @export
stat_num <- function (meta, min_nFeature_RNA = 200, max_nFeature_RNA = Inf,
                      percent.mt = 10, min.param.list = NULL, 
                      max.param.list = NULL,
                      str.param.list = NULL) 
{
    meta$Sample <- as.vector(meta$Sample)
    ave <- reshape2::melt(meta, id.var = "Sample")
    if (!is.null(str.param.list)) {
      ave <- ave[ -which(ave$variable %in% names(str.param.list)), ]
    }
    ave$value <- as.numeric(ave$value)
    stat <- ave %>% group_by(Sample, variable) %>% summarise(mean = mean(value))
    statd <- reshape2::dcast(stat, Sample ~ variable)
    statd$raw_percent.mt <- round(statd$percent.mt, 4)
    statd$raw_percent.ribo <- round(statd$percent.ribo, 4)
    statd$raw_nFeature <- round(statd$nFeature_RNA, 0)
    statd$raw_nCount <- round(statd$nCount_RNA, 0)
    statd$raw_number <- as.vector(table(meta$Sample))
    statd <- statd[, -which(colnames(statd) %in% c("percent.mt", 
        "percent.ribo", "nFeature_RNA", "nCount_RNA"))]
    Ave_OR_Total <- c("Ave_OR_Total", round(mean(statd$raw_percent.mt), 
        4), round(mean(statd$raw_percent.ribo), 4), round(mean(statd$raw_nFeature), 
        0), round(mean(statd$raw_nCount), 0), sum(statd$raw_number))
    statd <- rbind(statd, Ave_OR_Total)
    metaf <- meta[meta$nFeature_RNA > min_nFeature_RNA &
                    meta$percent.mt < percent.mt &
                    meta$nFeature_RNA < max_nFeature_RNA, ]
    fpos <- rep(TRUE, nrow(metaf))
    if (!is.null(min.param.list)) {
        minpos <- sapply(names(min.param.list), function(x) metaf[, 
            x] > min.param.list[[x]])
        minpos <- as.data.frame(minpos)
        fpos <- fpos & apply(minpos, 1, function(x) !sum(!x) > 0)
    }
    if (!is.null(max.param.list)) {
        maxpos <- sapply(names(max.param.list), function(x) metaf[, 
            x] < max.param.list[[x]])
        maxpos <- as.data.frame(maxpos)
        fpos <- fpos & apply(maxpos, 1, function(x) !sum(!x) > 0)
    }
    if (!is.null(str.param.list)) {
      strpos <- sapply(names(str.param.list), function(x) metaf[, 
                                                                x] == str.param.list[[x]])
      strpos <- as.data.frame(strpos)
      fpos <- fpos & apply(strpos, 1, function(x) !sum(!x) > 0)
      metaf <- metaf[, -which(colnames(metaf) %in% names(str.param.list))]
    }
    metaf <- metaf[fpos, ]
    ave <- reshape2::melt(metaf, id.var = "Sample")
    stat <- ave %>% group_by(Sample, variable) %>% summarise(mean = mean(value))
    statdf <- reshape2::dcast(stat, Sample ~ variable)
    statdf$filtered_percent.mt <- round(statdf$percent.mt, 4)
    statdf$filtered_percent.ribo <- round(statdf$percent.ribo, 4)
    statdf$filtered_nFeature <- round(statdf$nFeature_RNA, 0)
    statdf$filtered_nCount <- round(statdf$nCount_RNA, 0)
    statdf$filtered_number <- as.vector(table(metaf$Sample))
    statdf <- statdf[, -which(colnames(statdf) %in% c("percent.mt", 
        "percent.ribo", "nFeature_RNA", "nCount_RNA"))]
    Ave_OR_Total <- c("Ave_OR_Total", round(mean(statdf$filtered_percent.mt), 
        4), round(mean(statdf$filtered_percent.ribo), 4), round(mean(statdf$filtered_nFeature), 
        0), round(mean(statdf$filtered_nCount), 0), sum(statdf$filtered_number))
    statdf <- rbind(statdf, Ave_OR_Total)
    statmat <- left_join(statd, statdf, by = "Sample")
    statmat[is.na(statmat)] <- "-"
    return(statmat)
}

#' statistic of cell number by table().
#' 
#' @return data frame
#' 
#' @examples 
#' df <- statnum(meta200, rowname = 'Annotation', colname = 'Group')
#' 
#' @export
statnum <- function (metadf, rowname, colname, col_annotation = NULL) 
{
  tab <- table(metadf[, c(rowname, colname)])
  Total <- apply(tab, 1, sum)
  tab <- cbind(tab, Total)
  Total <- apply(tab, 2, sum)
  tab <- rbind(tab, Total)
  if (!is.null(col_annotation)) {
    unidf <- data.frame(unique(metadf[, c(colname, col_annotation)]), 
                        row.names = 1)
    tab <- rbind(cbind(t(unidf), "-"), tab)
    colnames(tab)[ncol(tab)] <- "Total"
  }
  return(tab)
}


#' Statistc of cell number before and after cell filtering.
#' 
#' @return A dataframe of cell number and other QC message.
#' 
#' @examples 
#' samples <- c(rep(c('sample1','sample2','sample3'),(ncol(pbmc) -1)/3),'sample1')
#' pbmc@meta.data$Sample <- samples
#' pbmc[["percent.ribo"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[SL]")
#' raw_sample <- pbmc$Sample
#' pbmc_filtered <- subset(pbmc, subset = nFeature_RNA > 800)
#' QC_table <- stat_table(pbmc_filtered, raw_sample)
#' 
#' @export
stat_table <- function (filtered_obj, raw_sample, percent.mt = "percent.mt", 
    percent.ribo = "percent.ribo", nFeature_RNA = "nFeature_RNA", 
    nCount_RNA = "nCount_RNA", Sample = "Sample") 
{
    meta <- filtered_obj@meta.data[, c(Sample, percent.mt, percent.ribo, 
        nFeature_RNA, nCount_RNA)]
    colnames(meta) <- c("Sample", "percent.mt", "percent.ribo", 
        "nFeature_RNA", "nCount_RNA")
    ave <- reshape2::melt(meta, id.var = "Sample")
    stat <- ave %>% group_by(Sample, variable) %>% summarise(mean = mean(value))
    statd <- dcast(stat, Sample ~ variable)
    statd$nFeature_RNA <- round(statd$nFeature_RNA, 0)
    statd$nCount_RNA <- round(statd$nCount_RNA, 0)
    statd$raw_number <- as.vector(table(raw_sample))
    statd$final_number <- as.vector(table(meta$Sample))
    statd$percent.mt <- round(statd$percent.mt, 4)
    statd$percent.ribo <- round(statd$percent.ribo, 4)
    Ave_OR_Total <- c("Ave_OR_Total", round(mean(statd$percent.mt), 
        4), round(mean(statd$percent.ribo), 4), round(mean(statd$nFeature_RNA), 
        0), round(mean(statd$nCount_RNA), 0), sum(statd$raw_number), 
        sum(statd$final_number))
    statdata <- rbind(statd, Ave_OR_Total)
    return(statdata)
}

#' QC statistic for one sample
#' 
#' @examples 
#' samples <- c(rep(c('sample1','sample2','sample3'),(ncol(pbmc) -1)/3),'sample1')
#' pbmc@@meta.data$Sample <- samples
#' pbmc[["percent.ribo"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[SL]")
#' raw_sample <- pbmc$Sample
#' pbmc_filtered <- subset(pbmc, subset = nFeature_RNA > 800)
#' QC_table_one <- stat_table_one_sample(pbmc_filtered, ncol(pbmc))
#' 
#' @export
stat_table_one_sample <- function (filtered_obj, raw_cell_num, percent.mt = "percent.mt", 
    percent.ribo = "percent.ribo", nFeature_RNA = "nFeature_RNA", 
    nCount_RNA = "nCount_RNA") 
{
    meta <- filtered_obj@meta.data[, c(percent.mt, percent.ribo, 
        nFeature_RNA, nCount_RNA)]
    colnames(meta) <- c("percent.mt", "percent.ribo", "nFeature_RNA", 
        "nCount_RNA")
    statd <- data.frame(nFeature_RNA = round(mean(meta$nFeature_RNA), 
        0), nCount_RNA = round(mean(meta$nCount_RNA), 0), percent.mt = round(mean(meta$percent.mt), 
        4), percent.ribo = round(mean(meta$percent.ribo), 4), 
        raw_number = raw_cell_num, final_number = ncol(filtered_obj), 
        stringsAsFactors = F)
    rownames(statd) <- c("Ave_OR_Total")
    return(statd)
}

#' Extract cells by proportion
#' 
#' @examples 
#' pbmc$Cellname <- colnames(pbmc)
#' subsetID(pbmc@@meta.data[,c('Annotation','Cellname')], num = 100)
#' 
#' @export
subsetID <- function (IDmat, num = NULL, fraction = NULL, min.cell = 10, 
                      seed = 123) 
{
  newI <- data.frame(ID = as.vector(IDmat[, 1]), Cellnames = as.vector(IDmat[,2]), 
                     stringsAsFactors = F)
  colnames(newI) <- c("ID", "Cellnames")
  mingroup <- min(table(newI$ID))
  if (is.null(num)) num <- min(table(IDmat$Annotation))
  frac <- round(num/nrow(newI), 5)
  set.seed(seed)
  if (frac > 1) {
    cell10000 <- newI
  }
  else if (mingroup * frac < 1) {
    cell10000 <- lapply(unique(newI$ID), function(x) {
      if (sum(newI$ID %in% x) < min.cell) {
        newI[newI$ID %in% x, ]
      }
      else if (sum(sample_frac(newI, frac)$ID %in% x) < 
               min.cell) {
        newI[newI$ID %in% x, ] %>% sample_n(min.cell)
      }
      else {
        newI[newI$ID %in% x, ] %>% sample_frac(frac)
      }
    })
    cell10000 <- do.call(rbind, cell10000)
  }
  else {
    cell10000 <- newI %>% group_by(ID) %>% sample_frac(frac)
  }
  if (!is.null(fraction)) {
    cell10000 <- newI %>% group_by(ID) %>% sample_frac(fraction)
    tab <- table(cell10000$ID)
    if (any(tab < min.cell) | length(unique(cell10000$ID)) != length(unique(newI$ID))) {
      min.group <- names(tab)[tab < min.cell]
      cell10000 <- cell10000[!cell10000$ID %in% min.group,]
      for ( i in min.group) {
        df <- newI[newI$ID == i,]
        set.seed(seed)
        if (nrow(df) > min.cell) df <- df[sample(1:nrow(df),min.cell),] 
        cell10000 <- rbind(cell10000, df)
      }
    }
  }
  return(cell10000)
}

#' Extract cells by expected number
#' 
#' @examples 
#' subsetID2(newI, expected.cell = 200, seed = 123)
#' 
#' @export
subsetID2 <- function (newI, expected.cell = 200, seed = 123) 
{
  colnames(newI) <- c("ID", "Cellnames")
  mingroup <- min(table(newI$ID))
  maxgroup <- max(table(newI$ID))
  if (mingroup > expected.cell) {
    set.seed(seed)
    cellext <- newI %>% group_by(ID) %>% sample_n(expected.cell)
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
    largegroupext <- largegroup %>% group_by(ID) %>% sample_n(expected.cell)
    cellext <- rbind(smallgroup, as.matrix(largegroupext))
  }
  return(cellext)
}


