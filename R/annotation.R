
#' Find candidate celltypes for each inputed group.
#'
#' @param object A seurat object.
#' @param active.assay Default assay in seurat object which used to model.
#' @param select.db Marker genes database for usage. Available options:
#' db_self(basic database), db_lung_zzm, db_immune_qss, db_cellmarker, panglaodb_hs,
#' db_ovarian_cancer_T;db_self_mm10(basic database), db_cellmarker_mm10, panglaodb_mm.
#' Default is basic database + lung_markers_zhangzemin_Cell2019.
#' @param add.db A dataframe. Add database by user. Colnames is c("cellName", "geneSymbol").
#' each geneset in "geneSymbol" is a character, such as 'CD3D,CD4,CCR7'.
#' @param DE The result of FindAllMarker as candidate genes.
#' @param do.DE Find defferent expression genes as candidate genes.
#' @param test.use Denotes which test to use if do.DE is TRUE.
#' Available options are:
#' \itemize{
#'  \item{"wilcox"} : Identifies differentially expressed genes between two
#'  groups of cells using a Wilcoxon Rank Sum test (default)
#'  \item{"t"} : Identify differentially expressed genes between two groups of
#'  cells using the Student's t-test.
#' }
#' @param core.num Number of cores used for do.DE. Default is 1
#' @param min.pct Only candidate genes that are detected in a minimum fraction of
#' min.pct cells. Default is 0.3
#' @param min.exp Only detected genes with expression level higher than min.exp.
#' @param show.top To show top n candidate predicted celltypes from the result of
#' GO enrichment. Default is 1
#' @param use.sigmoid.method Another available method for cadidate genes selection.
#' Using sigmoid function to identify candidate genes instead of min.pct.
#' @param sigmoid.min If use.sigmoid.method is TURE, only cells that are expressing
#' the detected gene in a minimum fraction of sigmoid.min sigmoid value.
#' @param sigmoid.pct If use.sigmoid.method is TURE, only genes that are candidate
#' genes in a minimum fraction of sigmoid.pct detected cells in one cluster.
#' @param HK.model.expr Another available method for cadidate genes selection. Get the
#' cutoff value base on the expression of house keep genes. The cutoff value from the
#' intersection point in a minimum fraction of area.frac area under curve
#' (dropOutRate x Expression).
#' @param HK.model.drop Another available method for cadidate genes selection. Get the
#' cutoff value base on the DropOutRate of house keep genes. The cutoff value from the
#' intersection point in a minimum fraction of area.frac area under curve
#' (dropOutRate x Expression).
#' @param HK.model.drop2 Another available method for cadidate genes selection.
#' The median value of the DropOutRate of house keep genes as cutoff value.
#' @param area.frac Minimum fraction of area under a curve from HK.model.
#' @param use.definite.markers Use strong marker to identify cell types. Default is F.
#' @param reduction Select the reduction method.
#' @param plot.result Plot predicted cell type and marker genes in GO terms.
#' @param species Optional speices: Homosapiens, Mouse.
#' @param mm10.to.human Transfer mouse homolog genes to human.
#' @param ... Arguments passed to signDE.
#' @param test.method Methods for hepergeometric test.
#' \itemize{
#'  \item{1}: lapply;
#'  \item{2}: phyper;
#'  \item{3}: Rcpp Ramanujan's approximation.
#'  }
#' @return A dataframe of celltype annotation.
#'
#' @examples
#' angrycell(pbmc) #pbmc is a seurat object.
#' angrycell(pbmc, min.pct = 0.3, show.top = 10, use.definite.markers = F,
#' reduction = 'umap')
#' angrycell(pbmc, min.pct = 0.3, show.top = 10, use.definite.markers = TRUE,
#' plot.result = TRUE, reduction = 'umap')
#' angrycell(pbmc, min.pct = 0.3, show.top = 10, use.sigmoid.method = TRUE,
#' sigmoid.min = 0.8, sigmoid.pct = 0.01)
#' angrycell(pbmc, HK.model.expr = TRUE, area.frac = 0.4)
#' angrycell(pbmc, HK.model.drop = TRUE, area.frac = 0.1)
#' angrycell(pbmc, HK.model.drop2 = TRUE)
#' @export
angrycell_v1.1.7 <- function(object = SeuratS4, select.db = c('db_self','db_lung_zzm'),
                             active.assay = "RNA", add.db = NULL, DE = NULL, do.DE = FALSE,
                             test.use = "wilcox",core.num = 1, min.pct = 0.3, min.exp = 0,
                             show.top = 1, use.sigmoid.method = FALSE,sigmoid.min = 0.85,
                             sigmoid.pct = 0.01, HK.model.expr = FALSE, HK.model.drop = FALSE,
                             HK.model.drop2 = FALSE, area.frac = 0.4,
                             use.definite.markers = FALSE, plot.result = FALSE,
                             reduction = 'tsne', species = 'Homosapiens', mm10.to.human = F,test.method = 2, ...) {
  object@active.assay <- active.assay
  mat <- object@assays[[active.assay]]@data
  # filter RP/MT genes
  mat <- mat[!grepl("^RP[SL]|^MT-", rownames(mat)), ]
  ident <- Idents(object)
  mat <- mat[, names(ident)]
  
  ### identify cell-type special geneset.
  if (!is.null(DE)) {
    genelist <- list()
    for (g in unique(DE$cluster)) {
      genelist[[g]] <- DE$gene[DE$cluster == g]
    }
  } else if (do.DE) {
    #use DE function in ./R/DE.R
    object@meta.data$Annotation <- paste0("cluster_", Idents(object))
    markerlist <- signDE(seurat = object, group.by = "Annotation",
                         BPPARAM = MulticoreParam(core.num), test.use = test.use,...)
    genelist <- list()
    for (g in names(markerlist)) {
      genelist[[g]] <- markerlist[[g]]$feature
    }
  } else if (use.sigmoid.method) {
    # Normalization optional mat <- apply(mat,2,function(x) (x-min(x)) /
    # (max(x) - min(x)))
    ##number of cells that sigmoid > 0.85 in each cluster.
    ##top_cluster: cluster x gene
    top_cluster <- apply(mat, 1, function(x) table(ident[sigmoid(x) >
                                                           sigmoid.min]))
    ## min.pct
    genelist <- list()
    cluster <- paste0("cluster", rownames(top_cluster))
    rownames(top_cluster) <- cluster
    genename <- colnames(top_cluster)
    identmat <- as.data.frame(table(ident))
    rownames(identmat) <- paste0("cluster", identmat$ident)
    identmat <- identmat[cluster, ]
    ##
    for (i in cluster) {
      genelist[[i]] <- genename[top_cluster[i, ] > sigmoid.pct *
                                  identmat[i, "Freq"]]
    }
  }
  ### 20200213 using housekeeeping genes and drop out rate to define
  ### expressing genes.
  else if (HK.model.expr) {
    genelist <- list()
    for (g in levels(ident)) {
      pos <- ident == g
      genem <- apply(mat[, pos], 1, mean)
      HKcutoff <- get_HKcutoff(object, active.assay, g, HKgene,
                               area.frac)[1]
      gpos <- genem > HKcutoff
      selected_gene <- rownames(mat)[gpos]
      cluster <- paste0("cluster", g)
      genelist[[cluster]] <- selected_gene
    }
  } else if (HK.model.drop) {
    genelist <- list()
    for (g in levels(ident)) {
      pos <- ident == g
      gened <- apply(mat[, pos], 1, function(x) sum(x != 0)/length(x))
      HKcutoff <- get_HKcutoff(object, active.assay, g, HKgene,
                               area.frac)[2]
      gpos <- gened > HKcutoff
      selected_gene <- rownames(mat)[gpos]
      cluster <- paste0("cluster", g)
      genelist[[cluster]] <- selected_gene
    }
  } else if (HK.model.drop2) {
    genelist <- list()
    for (g in levels(ident)) {
      pos <- ident == g
      gened <- apply(mat[, pos], 1, function(x) sum(x != 0)/length(x))
      # median drop out rate
      gpos <- gened > median(gened)
      selected_gene <- rownames(mat)[gpos]
      cluster <- paste0("cluster", g)
      genelist[[cluster]] <- selected_gene
    }
  } else {
    genelist <- list()
    for (g in levels(ident)) {
      pos <- ident == g
      gpos <- apply(mat[, pos], 1, function(x) sum(x > min.exp)/length(x) >
                      min.pct)
      selected_gene <- rownames(mat)[gpos]
      cluster <- paste0("cluster", g)
      genelist[[cluster]] <- selected_gene
    }
  }
  ## filtering common genes shared in all groups.
  comm <- Reduce(intersect, genelist)
  ##db
  if (species == 'Homosapiens') {
    Celltypes.db <- db_self
    if ( "db_lung_zzm" %in% select.db) {
      Celltypes.db <- rbind(db2, Celltypes.db)
    }
    if ( "db_immune_qss" %in% select.db) {
      Celltypes.db <- rbind(addm_qss, Celltypes.db)
    }
    if ( "db_cellmarker" %in% select.db) {
      Celltypes.db <- rbind(db1, Celltypes.db)
    }
    if ( "db_ovarian_cancer_T" %in% select.db) {
      Celltypes.db <- rbind(db3, Celltypes.db)
    }
  }else if (species == 'Mouse') {
    Celltypes.db <- db_self_mm10
    if ('db_cellmarker_mm10' %in% select.db) {
      Celltypes.db <- rbind(db4, Celltypes.db)
    }
  }else{
    stop('Only Homosapiens or Mouse!')
  }
  
  ## add annotation from user.
  if (!is.null(add.db)) {
    if (is.null(dim(add.db))) {
      add_type <- add.db
      names(add_type) <- colnames(Celltypes.db)
      Celltypes.db <- rbind(Celltypes.db, add_type)
    } else {
      add_type <- add.db
      colnames(add_type) <- colnames(Celltypes.db)
      Celltypes.db <- rbind(Celltypes.db, add_type)
    }
  }
  
  ## running celltypes identification
  N <- nrow(mat)
  celllist <- list()
  ## iterate each celltype
  for (n in names(genelist)) {
    candidate <- genelist[[n]]
    candidate <- candidate[!candidate %in% comm]
    if (mm10.to.human) {
      intergene <- intersect( homolog_Human_mouse[,2],candidate)
      candidate <- homolog_Human_mouse[intergene,1]
    }
    ## fisher's exact test
    if (test.method == 1) {
      termp <- get_terms_P1(candidate, N, Celltypes.db)
    }else if (test.method == 2) {
      termp <- get_terms_P2(candidate, N, Celltypes.db)
    }else if (test.method == 3) {
      termp <- get_terms_P3(candidate, N, Celltypes.db)
    }
    prdlist <- termp[1:show.top, c("Celltype", "Term.p.value", "BH.adj.p",
                                   "All_Term_Genes")]
    prdlist$Term.p.value <- round(as.double(prdlist$Term.p.value),
                                  5)
    prdlist$BH.adj.p <- round(as.double(prdlist$BH.adj.p), 5)
    Gene_in_Term <- c()
    termpct <- c()
    for (i in 1:show.top) {
      All_Term_Genes_set <- gsub(" ", "", strsplit(prdlist$All_Term_Genes[i],
                                                   ",")[[1]])
      Gene_in_Term_set <- intersect(candidate, All_Term_Genes_set)
      termpct <- c(termpct, paste0(length(Gene_in_Term_set), "/",
                                   length(All_Term_Genes_set)))
      Gene_in_Term_set <- paste(Gene_in_Term_set, collapse = ",")
      Gene_in_Term <- c(Gene_in_Term, Gene_in_Term_set)
    }
    orig_Idents <- gsub("cluster", "", n)
    prdlist <- data.frame(orig_Idents, prdlist, Gene_in_Term, termpct, stringsAsFactors = F)
    prdlist <- prdlist[, c(1, 2, 3, 4, 6, 7, 5)]
    rownames(prdlist) <- paste0("_", c(1:show.top))
    colnames(prdlist) <- c("Orig_Idents", "Celltype_predicted", "Term_Pvalue",
                           "BH_adj_Pvalue", "Gene_in_Term", "Term.pct", "All_Term_Genes")
    # use strong markers to identify celltype without runing GO enrichment
    if (use.definite.markers) {
      cellp <- list()
      for (ns in names(D_M)) {
        if (sum(D_M[[ns]] %in% candidate) == length(D_M[[ns]])) {
          cellp[[ns]] <- c(orig_Idents, ns, "-", "-", paste0(D_M[[ns]],collapse = ","),
                           1, paste0(D_M[[ns]], collapse = ","))
        }
      }
      for (ns in names(A_M)) {
        pos <- A_M[[ns]] %in% candidate
        comm_num <- sum(pos)
        if (comm_num > 0) {
          frac <- paste0(comm_num, "/", length(A_M[[ns]]))
          cellp[[ns]] <- c(orig_Idents, ns, "-", "-", paste0(A_M[[ns]][pos],collapse = ","),
                           frac, paste0(A_M[[ns]], collapse = ","))
        }
      }
      if (length(cellp) > 0) {
        cellpdata <- data.frame(t(data.frame(cellp)),stringsAsFactors = F)
        rownames(cellpdata) <- paste0("Strong_Marker", "_", 1:nrow(cellpdata))
        colnames(cellpdata) <- c("Orig_Idents", "Celltype_predicted",
                                 "Term_Pvalue", "BH_adj_Pvalue", "Gene_in_Term", "Term.pct",
                                 "All_Term_Genes")
        prdlist <- rbind(cellpdata, prdlist)
      }
    }
    celllist[[n]] <- prdlist
  }
  # merge annotation of each group.
  cellmat <- do.call(rbind, celllist)
  colnames(cellmat) <- c("Orig_Idents", "Celltype_predicted", "Term_Pvalue",
                         "BH_adj_Pvalue", "Gene_in_Term", "Term.pct", "All_Term_Genes")
  celldata <- as.data.frame(cellmat, stringsAsFactors = F)
  ## plot of cell and genes in GO terms
  if (plot.result) {
    show_result(object, active.assay, reduction, celldata, show.all = TRUE,
                return.grob = FALSE)
  }
  return(celldata)
}


#' Calculate p value of hypergeometric distribution test.
#'
#' @param candidate Genes to test.
#' @param N Number of cells.
#' @param db Marker genes database for usage.
#' @return Result of GO enrichment.
#' @export
get_terms_P1 <- function(candidate, N, db) {
  pvalue <- c()
  for (t in 1:nrow(db)) {
    term <- gsub(" ", "", strsplit(db$geneSymbol[t], ",")[[1]])
    n <- length(candidate)
    M <- length(term)
    k <- sum(candidate %in% term)
    p <- 1 - sum(sapply(0:k - 1, function(i) choose(n, i) * (M/N)^i *
                          (1 - M/N)^(n - i)))
    pvalue <- c(pvalue, p)
  }
  # add Benjamini–Hochberg procedure for Pvalue correction
  adjp <- p.adjust(pvalue, "BH")
  termp <- data.frame(as.vector(db$cellName), pvalue, adjp,
                      as.vector(db$geneSymbol), stringsAsFactors = F)
  termp <- termp[order(termp[, 2]), ]
  colnames(termp) <- c("Celltype", "Term.p.value", "BH.adj.p", "All_Term_Genes")
  return(termp)
}

#' Calculate p value of hypergeometric distribution test.
#' Hypergeometirc test by phyper function, which can get the same result as get_terms_P.
#'
#' @param candidate Genes to test.
#' @param N Number of cells.
#' @param db Marker genes database for usage.
#' @return Result of GO enrichment.
#' @export
get_terms_P2 <- function(candidate, N, db) {
  pvalue <- c()
  for (t in 1:nrow(db)) {
    term_gene <- gsub(" ", "", strsplit(db$geneSymbol[t], ",")[[1]])
    m = length(term_gene)
    n = N - m
    k = length(candidate)
    x = sum(candidate %in% term_gene)
    p <-  phyper(q=x - 1, m=m, n=n, k=k, lower.tail=FALSE)
    if (p < 0) {p = 0}
    pvalue <- c(pvalue, p)
  }
  adjp <- p.adjust(pvalue, "BH")
  termp <- data.frame(as.vector(db$cellName), pvalue, adjp,
                      as.vector(db$geneSymbol), stringsAsFactors = F)
  termp <- termp[order(termp[, 2]), ]
  colnames(termp) <- c("Celltype", "Term.p.value", "BH.adj.p", "All_Term_Genes")
  return(termp)
}


sigmoid <- function(x) {
  1/(1 + exp(-x))
}

#' Find the cutoff value of a curve.Default uses f(x) = ax^2 + bx +c to fit the curve.
#'
#' @param x x axis of the curve.
#' @param y y axis of the curve.
#' @param area.frac Fraction of area under the curve to find
#' @param ... other parameters in function trendline.
#' the intersection point on the curve.
#'
#' @return The upper and lower cutoff values
#' @export
get_cutoff <- function(x,y,area.frac,model='line3P',...){
  basicTrendline::trendline(x, y, model= model, ePos.x= 'top', linecolor
                            = 'red', CI.color = NA, xlab='PC', ylab='SD',...)
  tr <- basicTrendline::trendline_summary(x, y, model = model, summary = FALSE)
  a <- tr$parameter$a
  b <- tr$parameter$b
  c <- tr$parameter$c
  xt = -b/(2 * a)
  yt = (4 * a * c - b^2)/(4 * a)
  ## define x value of cross point if the peak is under x axis.
  if (yt < 0) {
    f3 <- function(x, a, b, c) a * x^2 + b * x + c
    result1 <- uniroot(f3, interval = c(0, xt), a = a, b = b, c = c,
                       tol = 1e-04)
    xt = result1$root
  }
  itg <- integrate(function(x) a * x^2 + b * x + c, 0, xt)
  area <- itg$value
  findprob <- function(f, interval, target) {
    optimize(function(x) {
      abs(integrate(f, 0, x)$value - target)
    }, interval)$minimum
  }
  mydensity <- function(x) a * x^2 + b * x + c
  expr_upper <- findprob(mydensity, interval = c(0, xt), target = area *
                           area.frac)
  
  findprob_dropoutRate <- function(f, interval, target) {
    optimize(function(x) {
      abs(integrate(f, 0, x)$value - x * (a * x^2 + b * x + c) -
            target)
    }, interval)$minimum
  }
  drop_upper_x <- findprob_dropoutRate(mydensity, interval = c(0, xt),
                                       target = area * area.frac)
  drop_upper_y <- a * drop_upper_x^2 + b * drop_upper_x + c
  c(expr_upper, drop_upper_y)
}
####nagative binomial model
##var and mean
varmean_fun <- function (x, phi) {
  return(x + phi*x^2)
}
##no-zero fraction and mean
prob_zero_fun <- function (mu, phi) {
  if (phi == 0) {
    return(exp(-mu))
  }
  phi_1 = 1 / phi
  return((phi_1 / (mu + phi_1)) ** phi_1)
}
#' get phi parameter of negative binomial model
#' @examples
#' get_phi(pbmc, 'RNA', 'Naive CD4 T')
#' @export
get_phi <- function(object, active.assay, celltype) {
  data <- GetAssayData(object, assay = active.assay, layer = 'counts')[,Idents(object) == celltype]
  Mean <-  as.vector(Matrix::rowMeans(data))
  Var <- as.vector(apply(data, 1, var))
  ##fit formula by python scipy.optimize.curve_fit
  scipy <- reticulate::import('scipy')
  fn <- reticulate::py_eval("lambda x,phi: x + phi*x**2")
  phi <- scipy$optimize$curve_fit(fn,Mean,Var)
  return(phi[[1]])
}
#' get zero-fraction cut off value from negative binomial model.
#' Cluster with large phi is supposed to contain multi-cell types,
#' and cell marker genes have high zero fraction.
#' @examples
#' get_cutoff_NB(pbmc, 'RNA', 'CD8 T',0.3)
#' @export
get_cutoff_NB <- function(object, active.assay, celltype, initial.zero.frac) {
  phi <- get_phi(object, active.assay, celltype)
  if (phi < 0.01) {
    drop.cut <- 0.6
    return(drop.cut)
  }else if (phi > 1) {
    drop.cut <- 0.2
    return(drop.cut)
  }
  scalel <- log10(10) - log10(0.01)
  drop.cut <- initial.zero.frac * (log10(10) - log10(phi))/ scalel * 2
  return(as.numeric(drop.cut))
}
#' gene list for each cluster
#' @examples
#'
#' @export
get_genelist <- function(object, active.assay, celltype, min.pct = 0.3) {
  data <- object@assays[[active.assay]]@data[,Idents(object) == celltype]
  dropOutRate <- apply(data, 1, function(x) {
    sum(x == 0)/length(x)
  })
  if (min.pct == 'auto') {
    min.pct <- get_cutoff_NB(object, active.assay, celltype, initial.zero.frac = 0.3)
  }
  genelist <- rownames(data)[1 - dropOutRate > min.pct]
  return(genelist)
}
#' Using housekeeeping genes and negative binomial model to get the cutoff value
#' of dropout rate.
#'
#' @param object A seurat object.
#' @param active.assay Default assay in seurat object which used to model.
#' @param celltype The cell cluster name.
#' @param HKgenes House Keeping genes to model.
#' @param min.pct cell Fraction of no-zero expression.
#' the intersection point on the curve.
#'
#' @return gene list
#' @export
get_HKcutoff <- function(object, active.assay, celltype, HKgenes, min.pct) {
  HKgene_indata <- rownames(object)[rownames(object) %in% HKgenes]
  subs <- object[HKgene_indata,]
  data <- subs@assays[[active.assay]]@data
  # baseg <- rownames(data)[rownames(data) %in% HKgenes]
  datacl <- data[,Idents(object) == celltype]
  dropOutRate <- apply(datacl, 1, function(x) {
    sum(x == 0)/length(x)
  })
  meanExp <- apply(datacl, 1, mean)
  ##HKgene predicted drop out rate from NB model
  phi <- get_phi(subs, active.assay, celltype)
  pred_drop <- prob_zero_fun(meanExp, phi)
  ##filtering gene, to get genes above the HKgene NB curve.
  gene_pos <- dropOutRate > min.pct & dropOutRate > pred_drop
  geneset <- rownames(datacl)[gene_pos]
  return(geneset)
}



#' Find candidate celltypes for each inputed group.
#'
#' @param object A seurat object.
#' @param active.assay Default assay in seurat object which used to model.
#' @param select.db Marker genes database for usage. Available options:
#' db_self(basic database), db_lung_zzm, db_immune_qss, db_cellmarker,
#' db_ovarian_cancer_T;db_self_mm10(basic database), db_cellmarker_mm10.
#' Default is basic database + lung_markers_zhangzemin_Cell2019.
#' @param add.db A dataframe. Add database by user. Colnames is c("cellName", "geneSymbol").
#' each geneset in "geneSymbol" is a character, such as 'CD3D,CD4,CCR7'.
#' @param layers select different scale of annotation.
#' e.t. 1: Lymphoid; 2: T cell; 3: CD4+ naive T cell.
#' @param validate.cell To filter candidate cell types by expression of specific markers.
#' @param DE The result of FindAllMarker as candidate genes.
#' @param do.DE Find defferent expression genes as candidate genes.
#' @param test.use Denotes which test to use if do.DE is TRUE.
#' Available options are:
#' \itemize{
#'  \item{"wilcox"} : Identifies differentially expressed genes between two
#'  groups of cells using a Wilcoxon Rank Sum test (default)
#'  \item{"t"} : Identify differentially expressed genes between two groups of
#'  cells using the Student's t-test.
#' }
#' @param core.num Number of cores used for do.DE. Default is 1
#' @param min.pct Only candidate genes that are detected in a minimum fraction of
#' min.pct cells. 0.3 is recommended. Default is auto.
#' @param min.exp Only detected genes with expression level higher than min.exp.
#' @param show.top To show top n candidate predicted celltypes from the result of
#' GO enrichment. Default is 1
#' @param use.sigmoid.method Another available method for cadidate genes selection.
#' Using sigmoid function to identify candidate genes instead of min.pct.
#' @param sigmoid.min If use.sigmoid.method is TURE, only cells that are expressing
#' the detected gene in a minimum fraction of sigmoid.min sigmoid value.
#' @param sigmoid.pct If use.sigmoid.method is TURE, only genes that are candidate
#' genes in a minimum fraction of sigmoid.pct detected cells in one cluster.
#' @param HK.model.drop Another available method for cadidate genes selection. Get the
#' cutoff value base on the DropOutRate of house keep genes. The cutoff value from the
#' intersection point in a minimum fraction of area.frac area under curve
#' (dropOutRate x Expression).
#' @param HK.model.drop2 Another available method for cadidate genes selection.
#' The median value of the DropOutRate of house keep genes as cutoff value.
#' @param use.definite.markers Use strong marker to identify cell types. Default is F.
#' @param reduction Select the reduction method.
#' @param plot.result Plot predicted cell type and marker genes in GO terms.
#' @param species Optional speices: Homosapiens, Mouse.
#' @param mm10.to.human Transfer mouse homolog genes to human.
#' @param test.method Methods for hepergeometric test.
#' \itemize{
#'  \item{1}: lapply;
#'  \item{2}: phyper;
#'  \item{3}: Rcpp Ramanujan's approximation.
#'  }
#' @param ... Arguments passed to anno_openai.
#' 
#' @return A dataframe of celltype annotation.
#'
#' @examples
#' angrycell(pbmc) #pbmc is a seurat object.
#' angrycell(pbmc, min.pct = 0.3, show.top = 10, layer = c(1),
#' reduction = 'umap') #Annotated by main groups.
#' angrycell(pbmc, min.pct = 0.3, show.top = 10,layer = c(2,3), validate.cell = F,
#' reduction = 'umap') # no validation of predicted cell types.
#' angrycell(pbmc, min.pct = 0.3, show.top = 10,layer = c(2,3), validate.cell = T,
#' plot.result = TRUE, reduction = 'umap')
#' angrycell(pbmc, min.pct = 0.3, show.top = 10, use.sigmoid.method = TRUE,
#' sigmoid.min = 0.8, sigmoid.pct = 0.01)
#' angrycell(pbmc, HK.model.expr = TRUE, area.frac = 0.4)
#' angrycell(pbmc, HK.model.drop = TRUE, area.frac = 0.1)
#' angrycell(pbmc, HK.model.drop2 = TRUE)
#' 
#' @export
angrycell <- function(object, select.db = c('db_self'),
                      layers = c(2,3), validate.cell = F,
                      active.assay = "RNA", add.db = NULL, DE = NULL, do.DE = FALSE,
                      test.use = "wilcox",core.num = 1, min.pct = 'auto', min.exp = 0,
                      show.top = 1, use.sigmoid.method = FALSE,sigmoid.min = 0.85,
                      sigmoid.pct = 0.01, HK.model.drop = FALSE,
                      HK.model.drop2 = FALSE, return.genelist = F,
                      use.definite.markers = FALSE, plot.result = FALSE,
                      reduction = 'tsne', species = 'Homosapiens', mm10.to.human = F,
                      test.method = 2, ...) {
  options(stringsAsFactors = F)
  object@active.assay <- active.assay
  mat <- GetAssayData(object, assay = active.assay)
  # filter RP/MT genes
  mat <- mat[!grepl("^RP[SL]|^MT-", rownames(mat)), ]
  ident <- Idents(object)
  mat <- mat[, names(ident)]
  #initial min.pct
  if (min.pct == 'auto') {
    min.pct.list <- sapply(levels(ident), function(x) {
      get_cutoff_NB(object, active.assay, x, initial.zero.frac = 0.3)
    })
    names(min.pct.list) <- levels(ident)
  }else{
    min.pct.list <- rep(min.pct, length(levels(ident)))
    names(min.pct.list) <- levels(ident)
  }
  ### identify cell-type special geneset.
  if (!is.null(DE)) {
    genelist <- list()
    for (g in unique(DE$cluster)) {
      genelist[[g]] <- DE$gene[DE$cluster == g]
    }
  } else if (do.DE) {
    #use DE function in ./R/DE.R
    object@meta.data$Annotation <- paste0("cluster_", Idents(object))
    markerlist <- signDE(seurat = object, group.by = "Annotation",
                         BPPARAM = MulticoreParam(core.num), test.use = test.use)
    genelist <- list()
    for (g in names(markerlist)) {
      genelist[[g]] <- markerlist[[g]]$feature
    }
  } else if (use.sigmoid.method) {
    # Normalization optional mat <- apply(mat,2,function(x) (x-min(x)) /
    # (max(x) - min(x)))
    ##number of cells that sigmoid > 0.85 in each cluster.
    ##top_cluster: cluster x gene
    top_cluster <- apply(mat, 1, function(x) table(ident[sigmoid(x) >
                                                           sigmoid.min]))
    ## min.pct
    genelist <- list()
    cluster <- paste0("cluster", rownames(top_cluster))
    rownames(top_cluster) <- cluster
    genename <- colnames(top_cluster)
    identmat <- as.data.frame(table(ident))
    rownames(identmat) <- paste0("cluster", identmat$ident)
    identmat <- identmat[cluster, ]
    ##
    for (i in cluster) {
      genelist[[i]] <- genename[top_cluster[i, ] > sigmoid.pct *
                                  identmat[i, "Freq"]]
    }
  }
  ### 20200213 using housekeeeping genes and drop out rate to define
  ### expressing genes.
  else if (HK.model.drop) {
    genelist <- list()
    for (g in levels(ident)) {
      min.pct <- min.pct.list[names(min.pct.list) == g]
      #get the gene set by HKgene NB model
      HKcutoff <- get_HKcutoff(object, active.assay, g, HKgene,
                               min.pct)
      cluster <- paste0("cluster", g)
      genelist[[cluster]] <- HKcutoff
    }
  } else if (HK.model.drop2) {
    genelist <- list()
    for (g in levels(ident)) {
      pos <- ident == g
      gened <- apply(mat[, pos], 1, function(x) sum(x != 0)/length(x))
      # median drop out rate
      gpos <- gened > median(gened)
      selected_gene <- rownames(mat)[gpos]
      cluster <- paste0("cluster", g)
      genelist[[cluster]] <- selected_gene
    }
  } else {
    genelist <- list()
    for (g in levels(ident)) {
      pos <- ident == g
      min.pct <- min.pct.list[names(min.pct.list) == g]
      gpos <- apply(mat[, pos], 1, function(x) sum(x > min.exp)/length(x) >
                      min.pct)
      selected_gene <- rownames(mat)[gpos]
      cluster <- paste0("cluster", g)
      genelist[[cluster]] <- selected_gene
    }
    if (return.genelist) {
      return(genelist)
    }
  }
  ### 20250120 add openai function
  if(select.db == 'openai') {
    annodf <- anno_ellmer(genelist = genelist, llm_function = 'openai', ...)
    return(annodf)
  }
  ### 20250120 End.
  ## filtering common genes shared in all groups.
  comm <- Reduce(intersect, genelist)
  ##db
  if (species == 'Homosapiens') {
    Celltypes.db <- db_self
    if (select.db == 'db_self') {
      if (sum(!layers %in% c(1,2,3)) > 0 ) {
        stop('You must only select cell layers from 1,2,3. For example, layers = c(1,2).')
      }
      Celltypes.db <- do.call(rbind,dblayers[layers])
    }
    if ( "db_lung_zzm" %in% select.db) {
      Celltypes.db <- rbind(db2, Celltypes.db)
    }
    if ( "db_immune_qss" %in% select.db) {
      Celltypes.db <- rbind(addm_qss, Celltypes.db)
    }
    if ( "db_cellmarker" %in% select.db) {
      Celltypes.db <- rbind(db1, Celltypes.db)
    }
    if ( "db_ovarian_cancer_T" %in% select.db) {
      Celltypes.db <- rbind(db3, Celltypes.db)
    }
  }else if (species == 'Mouse') {
    Celltypes.db <- db_self_mm10
    if ('db_cellmarker_mm10' %in% select.db) {
      Celltypes.db <- rbind(db4, Celltypes.db)
    }
  }else{
    stop('Only Homosapiens or Mouse!')
  }
  
  ## add annotation from user.
  if (!is.null(add.db)) {
    if (is.null(dim(add.db))) {
      add_type <- add.db
      names(add_type) <- colnames(Celltypes.db)
      Celltypes.db <- rbind(Celltypes.db, add_type)
    } else {
      add_type <- add.db
      colnames(add_type) <- colnames(Celltypes.db)
      Celltypes.db <- rbind(Celltypes.db, add_type)
    }
  }
  
  ## running celltypes identification
  N <- nrow(mat)
  celllist <- list()
  ## iterate each celltype
  for (n in names(genelist)) {
    candidate <- genelist[[n]]
    candidate <- candidate[!candidate %in% comm]
    if (mm10.to.human) {
      intergene <- intersect(homolog_Human_mouse[,2],candidate)
      candidate <- homolog_Human_mouse[intergene,1]
    }
    ## fisher's exact test
    if (test.method == 1) {
      termp <- get_terms_P1(candidate, N, Celltypes.db)
    }else if (test.method == 2) {
      termp <- get_terms_P2(candidate, N, Celltypes.db)
    }
    ### consider layers or not
    if (select.db == 'db_self') {
      ##filter cell types by specific markers and add reference
      specificmat <- dbspf[termp$Celltype,]
      refmat <- dbref[termp$Celltype,]
      prdlist <- cbind(termp,
                       specificmat[,'Specific_marker', drop = F],
                       refmat[,'Reference', drop = F])
      if (validate.cell) {
        pos <- specificmat$Specific_marker %in% candidate
        if (sum(pos) == 0) {
          warning(paste0('There is not specific marker in the expressed marker list of ',n,
                         '. Return the first line for default.'))
          pos = 1
        }
        prdlist <- prdlist[pos,]
      }
      ranges <- min(nrow(prdlist),show.top)
      prdlist <- prdlist[1:ranges, c("Celltype", "Term.p.value", "BH.adj.p",
                                     "All_Term_Genes",'Specific_marker','Reference')]
      prdlist$Term.p.value <- round(as.double(prdlist$Term.p.value),
                                    5)
      prdlist$BH.adj.p <- round(as.double(prdlist$BH.adj.p), 5)
      Gene_in_Term <- c()
      termpct <- c()
      for (i in 1:ranges) {
        All_Term_Genes_set <- gsub(" ", "", strsplit(prdlist$All_Term_Genes[i],
                                                     ",")[[1]])
        Gene_in_Term_set <- intersect(candidate, All_Term_Genes_set)
        termpct <- c(termpct, paste0(length(Gene_in_Term_set), "/",
                                     length(All_Term_Genes_set)))
        Gene_in_Term_set <- paste(Gene_in_Term_set, collapse = ",")
        Gene_in_Term <- c(Gene_in_Term, Gene_in_Term_set)
      }
      orig_Idents <- gsub("cluster", "", n)
      prdlist <- data.frame(orig_Idents, prdlist, Gene_in_Term, termpct, stringsAsFactors = F)
      prdlist <- prdlist[, c(1, 2, 6, 3, 4, 8, 9, 5, 7)]
      rownames(prdlist) <- paste0("_", c(1:ranges))
      colnames(prdlist) <- c("Orig_Idents", "Celltype_predicted", 'Specific_marker', "Term_Pvalue",
                             "BH_adj_Pvalue", "Gene_in_Term", "Term.pct", "All_Term_Genes",
                             'Reference')
    }else{
      
      ## no layer method
      prdlist <- termp[1:show.top, c("Celltype", "Term.p.value", "BH.adj.p",
                                     "All_Term_Genes")]
      prdlist$Term.p.value <- round(as.double(prdlist$Term.p.value),
                                    5)
      prdlist$BH.adj.p <- round(as.double(prdlist$BH.adj.p), 5)
      Gene_in_Term <- c()
      termpct <- c()
      for (i in 1:show.top) {
        All_Term_Genes_set <- gsub(" ", "", strsplit(prdlist$All_Term_Genes[i],
                                                     ",")[[1]])
        Gene_in_Term_set <- intersect(candidate, All_Term_Genes_set)
        termpct <- c(termpct, paste0(length(Gene_in_Term_set), "/",
                                     length(All_Term_Genes_set)))
        Gene_in_Term_set <- paste(Gene_in_Term_set, collapse = ",")
        Gene_in_Term <- c(Gene_in_Term, Gene_in_Term_set)
      }
      orig_Idents <- gsub("cluster", "", n)
      prdlist <- data.frame(orig_Idents, prdlist, Gene_in_Term, termpct, stringsAsFactors = F)
      prdlist <- prdlist[, c(1, 2, 3, 4, 6, 7, 5)]
      rownames(prdlist) <- paste0("_", c(1:show.top))
      colnames(prdlist) <- c("Orig_Idents", "Celltype_predicted", "Term_Pvalue",
                             "BH_adj_Pvalue", "Gene_in_Term", "Term.pct", "All_Term_Genes")
      # use strong markers to identify celltype without runing GO enrichment
      if (use.definite.markers) {
        cellp <- list()
        for (ns in names(D_M)) {
          if (sum(D_M[[ns]] %in% candidate) == length(D_M[[ns]])) {
            cellp[[ns]] <- c(orig_Idents, ns, "-", "-", paste0(D_M[[ns]],collapse = ","),
                             1, paste0(D_M[[ns]], collapse = ","))
          }
        }
        for (ns in names(A_M)) {
          pos <- A_M[[ns]] %in% candidate
          comm_num <- sum(pos)
          if (comm_num > 0) {
            frac <- paste0(comm_num, "/", length(A_M[[ns]]))
            cellp[[ns]] <- c(orig_Idents, ns, "-", "-", paste0(A_M[[ns]][pos],collapse = ","),
                             frac, paste0(A_M[[ns]], collapse = ","))
          }
        }
        if (length(cellp) > 0) {
          cellpdata <- data.frame(t(data.frame(cellp)),stringsAsFactors = F)
          rownames(cellpdata) <- paste0("Strong_Marker", "_", 1:nrow(cellpdata))
          colnames(cellpdata) <- c("Orig_Idents", "Celltype_predicted",
                                   "Term_Pvalue", "BH_adj_Pvalue", "Gene_in_Term", "Term.pct",
                                   "All_Term_Genes")
          prdlist <- rbind(cellpdata, prdlist)
        }
      }
    }
    
    celllist[[n]] <- prdlist
  }
  # merge annotation of each group.
  cellmat <- do.call(rbind, celllist)
  # colnames(cellmat) <- c("Orig_Idents", "Celltype_predicted", "Term_Pvalue",
  #                       "BH_adj_Pvalue", "Gene_in_Term", "Term.pct", "All_Term_Genes")
  celldata <- as.data.frame(cellmat, stringsAsFactors = F)
  ## plot of cell and genes in GO terms
  if (plot.result) {
    show_result(object, active.assay, reduction, celldata, show.all = TRUE,
                return.grob = FALSE)
  }
  return(celldata)
}

################## Annotation
#' Annotate raw clusters by max expression, and return an object with annotation Idents
#' @examples
#' markerdf <- data.frame(celltypes = c('T','T','NK','NK','Mono', 'Mono','DC','DC','B','B','Platelet'), 
#' markers = c('CD3D','CD3E','NCAM1','NKG7','CD14','FCGR3A','CST3','CD1C','CD79A','MS4A1','PPBP'))
#' colnames(markerdf)
#' #[1] "celltypes" "markers"
#' pbmc_anno <- anno_top_gene(pbmc, markerdf, group.by = 'seurat_clusters')
#' DimPlot(pbmc_anno, label = T)
#' 
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
#' 
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
#' 
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
#' 
#' @export
anno_ellmer <- function(deg = NULL, genelist = NULL, tissuename = NULL,
                        base_url = "http://chatapi.littlewheat.com/v1",
                        api_key = 'sk-HgtySiUAhSLiZTlDRhNE7aEbERJOuSumUveDxYfAUy8YvDfM',
                        model = "gpt-3.5-turbo", 
                        llm_function = c('openai','ollama'),
                        ollama_model = 'llama3.2',
                        return.content = F,
                        seed = 1234) {
  llm_function <- match.arg(NULL, choices = llm_function)
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
  content = paste(c(list(paste0("Identify cell types of ", 
                         tissuename, " cells using the following markers separately for each row."),
                  " Only provide the cell type name.",
                  " Some can be a mixture of multiple cell types."),
                  input), collapse = "\n")
  if (return.content) return(content)
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
      content <- paste0("Identify cell types of ", 
                        tissuename, " cells using the following markers separately for each row.",
      " Only provide the cell type name. Do not interpret.",
      " Some can be a mixture of multiple cell types.")
      for (i in input) {
        text <- chat$chat(paste0(content))
      }
  }
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
#' 
#' @export
annocell <- function(object, species, assay = 'RNA', raw_cluster = NULL, 
                     method = c('SingleR','AUCell','topgene','angrycell'), 
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
                     min.pct = 0.3,
                     db = NULL,
                     n.cores = 4,
                     ...) {
  method <- match.arg(NULL, choices = method)
  ref.label <-  match.arg(NULL, choices = ref.label)
  scsig.subset <-  match.arg(NULL, choices = scsig.subset)
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
  if (method == 'topgene') {
    object <- anno_top_gene(object, markerdf, group.by = raw_cluster)
    object@meta.data[,label_raw_cluster_colname] <- Idents(object)
  }
  if (method == 'angrycell') {
    if (is.null(db)) db <- ifelse(species == 'Human', 'panglaodb_hs', 'panglaodb_mm')
    if (!is.null(raw_cluster)) Idents(object) <- object@meta.data[, raw_cluster]
    df <- angrycell(object, select.db = db, show.top = 1, min.pct = min.pct, core.num = n.cores,...)
    new.id <- df$Celltype_predicted
    names(new.id) <- df$Orig_Idents
    Idents(object) <- object@meta.data[,raw_cluster] 
    object <- RenameIdents(object, new.id)
    object@meta.data[,label_raw_cluster_colname] <- Idents(object)
  }
  return(object)
}



