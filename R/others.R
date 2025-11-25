#' In order to avoid memory problem, Convert large dgCMatrix to Matrix by separation.
#'
#' @export
bigdataTOmat <- function(bigdata, min_num = 10000, transpose = F,
                         as.type = c('matrix', 'TsparseMatrix', 'dgCMatrix')) {
  as.type <- match.arg(arg = as.type, choices = c('matrix', 'TsparseMatrix', 'dgCMatrix'))
  matlist <- list()
  colnum <- ncol(bigdata)
  if (colnum < min_num){
    if (transpose) {
      matlist[['mat1']] <- as(Matrix::t(bigdata), as.type)
    }else{
      matlist[['mat1']] <- as(bigdata, as.type)
    }
    
  }else{
    message(paste0('The Matrix is too large. Divide the Matrix into smaller parts
            with ',min_num,' cells and convert them to matrices respectively.'))
    mod <- colnum %/% min_num
    matlist <- pbapply::pblapply(1:mod, function(x) {
      mat <- bigdata[,(min_num*(x-1) + 1):(min_num*x)]
      if (transpose) mat <- Matrix::t(mat)
      mat <- as(mat, as.type)
      return(mat)
    }
    )
    residue <- colnum %% min_num
    if (residue != 0) {
      mat <- bigdata[,(colnum - residue + 1) : colnum]
      if (transpose) mat <- Matrix::t(mat)
      matlist[['last_mat']] <- as(mat, as.type)
    }
  }
  if (as.type == 'matrix') matrix_mat <- do.call(cbind, matlist)
  if (as.type == 'TsparseMatrix') matrix_mat <- Matrix::.bdiag(matlist)
  if (as.type == 'dgCMatrix') matrix_mat <- Matrix::bdiag(matlist)
  return(matrix_mat)
}



mark_order_to_gene <- function (namelist) {
  no_numeric_prefix <- gsub("[0-9]+$", "", namelist)
  numeric_suffix <- suppressWarnings(as.numeric(gsub("(^.*[[:alpha:]])([0-9]+)$", 
                                                     "\\2", namelist)))
  numeric_suffix[is.na(numeric_suffix)] <- 1
  order <- sprintf("%08d", numeric_suffix)
  tempgenename <- paste0(no_numeric_prefix, "_", order, "_order_", 
                         namelist)
  return(tempgenename)
}
#' Order genes for plot.
#' 
#' @param method DEG or seriation. DEG: Order genes by DEG; seriation: Order genes by R package seriation. Default is DEG.
#' @param gene.Annot Data.frame, included 2 columns: annotation and gene. If sort.name.per.group.byAnnot is true, order genes by gene.Annot data frame in each cluster.
#' @param seriate.method Methods only for seriation. Options: BEA_TSP, BEA, PCA, PCA_angle. Default is BEA_TSP.
#' 
#' @examples 
#' marker_to_plot <- c('CD3D','CD4','CD8A','CCR7','IL7R','NCAM1','CD1C','CD68','CD14','CD79A','PPBP','FCGR3A')
#' ordered_gene <- order_heatmap_gene(pbmc, marker_to_plot, method = 'DEG')
#' DotPlot(pbmc, features = ordered_gene)
#' 
#' @export
order_heatmap_gene <- function (object, features, assay = "RNA", slot = "data", min.pct = 0, 
    method = "DEG", sort.name.per.group = F, name.decreasing = F, 
    sort.name.per.group.byAnnot = F, gene.Annot = NULL, seriate.method = c("BEA_TSP", 
        "BEA", "PCA", "PCA_angle"), ident = Idents(object), cell.order = NULL) 
{
    seriate.method <- base::match.arg(arg = seriate.method, choices = c("BEA_TSP", 
        "BEA", "PCA", "PCA_angle"), several.ok = T)
    Idents(object) <- ident
    if (!is.null(cell.order)) {
        Idents(object) <- factor(ident, levels = cell.order)
    }
    pctmat <- apply(GetAssayData(object, assay = assay, layer = slot)[features, 
        ], 1, function(x) {
        exp <- data.frame(cluster = ident, expression = x)
        sapply(levels(ident), function(x) {
            series <- exp[exp[, "cluster"] == x, "expression"]
            sum(series > 0)/length(series)
        })
    })
    pct <- apply(pctmat, 2, max)
    if (any(pct < min.pct)) {
        message(paste0("Genes with low percentage of expression were removed: \n", 
            paste(features[pct < min.pct], collapse = ",")))
        features <- features[pct >= min.pct]
    }
    if (method == "DEG") {
        data <- GetAssayData(object, assay = assay, layer = slot)[features, ]
        subs <- CreateSeuratObject(counts = data, meta.data = object@meta.data)
        subs@assays[[assay]]$data <- data
        Idents(subs) <- Idents(object)
        sDEG <- FindAllMarkers(subs, assay = assay, slot = slot, 
            only.pos = T, logfc.threshold = 0, min.pct = min.pct)
        pos <- !features %in% unique(sDEG$gene)
        if (any(pos)) {
            warning(paste("Gene not in DEGlist: ", paste(features[pos], 
                collapse = ",")))
        }
        if (sort.name.per.group) {
            temp_gene_name <- mark_order_to_gene(sDEG$gene)
            sDEG$gene <- temp_gene_name
            if (name.decreasing) {
                sDEG <- sDEG %>% dplyr::group_by(cluster) %>% 
                  arrange(desc(gene), .by_group = T)
            }
            else {
                sDEG <- sDEG %>% dplyr::group_by(cluster) %>% 
                  arrange(gene, .by_group = T)
            }
            sDEG$gene <- gsub("^.*_order_", "", sDEG$gene)
        }
        if (sort.name.per.group.byAnnot) {
            colnames(gene.Annot) <- c("annotation", "gene")
            sDEG <- left_join(sDEG, gene.Annot, by = "gene")
            dupgene <- sDEG$gene[duplicated(sDEG$gene)]
            df <- sDEG[, c("cluster", "gene")]
            sDEG <- sDEG %>% dplyr::group_by(cluster) %>% arrange(desc(annotation), 
                .by_group = T)
            sDEG <- sDEG[!duplicated(sDEG$gene), ]
            sDEG$annotation <- NULL
        }
        topgene <- sDEG
        topgene$index <- 1:nrow(topgene)
	pos <- colnames(topgene) %in% "avg_log2FC"
	if (any(pos)) {colnames(topgene)[pos] <- 'avg_logFC'}
        mat <- topgene[, c("index", "gene", "avg_logFC", "pct.1")]
        mat <- mat[order(mat$gene, mat$avg_logFC, mat$pct.1, 
            decreasing = T), ]
        mat <- mat[!duplicated(mat$gene), ]
        index <- sort(mat$index)
        selected_genes <- topgene$gene[index]
        return(selected_genes)
    }
    else if (method == "seriation") {
        orderlist <- list()
        aveobjs <- AverageExpression(object = object, assays = assay, 
            layer = slot, features = features, return.seurat = T)
        mat <- GetAssayData(aveobjs, assay = assay, slot = slot)
        mat <- as.matrix(mat)
        if (min(mat) < 0) {
            mat <- max(mat) - mat
        }
        o = seriation::seriate(mat, method = seriate.method)
        orderlist[["gene_order"]] <- seriation::get_order(o, 1)
        orderlist[["cell_order"]] <- seriation::get_order(o, 2)
        return(orderlist[["gene_order"]])
    }
}

#' adjust length of annotation text in heatmap.
#' 
#' @export
adj_string_len <- function (x, num) 
{
  out <- c()
  for (i in 1:length(x)) {
    times <- nchar(x[i])%/%num
    if (times != 0) {
      a <- strsplit(x[i], "")[[1]]
      index = 0
      for (j in 1:times) {
        index <- index + num
        a <- append(a, "\n", after = index)
        index + 1
      }
      a <- paste0(a, collapse = "")
    }
    else {
      a <- x[i]
    }
    out <- c(out, a)
  }
  return(out)
}

#' Unknown function
FoldChange <- function (data1, data2, pseudocount.use = 1, base = exp(1)) {
  log_ave1 <- log(x = rowMeans(x = expm1(x = data1)) + pseudocount.use, 
                  base = base)
  log_ave2 <- log(x = rowMeans(x = expm1(x = data2)) + pseudocount.use, 
                  base = base)
  fc <- log_ave1 - log_ave2
  return(fc)
}

#' Unknown function
FoldChange.Assay <- function (object, annotation, ident.1, ident.2,
                              features = NULL, 
          slot = "data", pseudocount.use = 1, fc.name = NULL, base = 2, 
          ...) 
{
  data <- GetAssayData(object = object, slot = slot)
  mean.fxn <- switch(EXPR = slot, data = function(x) {
    return(log(x = rowMeans(x = expm1(x = x)) + pseudocount.use, 
               base = base))
  }, scale.data = rowMeans, function(x) {
    return(log(x = rowMeans(x = x) + pseudocount.use, base = base))
  })
  base.text <- ifelse(test = base == exp(1), yes = "", no = base)
  if (is.null(fc.name)) {
    fc.name <- ifelse(test = slot == "scale.data", yes = "avg_diff", 
                      no = paste0("avg_log", base.text, "FC"))
  }
  cells.1 <- colnames(object)[object@meta.data[, annotation] %in% 
                                ident.1]
  cells.2 <- colnames(object)[object@meta.data[, annotation] %in% 
                                ident.2]
  FoldChange.default(object = data, cells.1 = cells.1, cells.2 = cells.2, 
                     features = features, mean.fxn = mean.fxn, fc.name = fc.name)
}

#' Unknown function
FoldChange.default <- function (object, cells.1, cells.2, mean.fxn, fc.name, features = NULL) 
{
  if (is.null(features)) {
    features <- rownames(x = object)
  }
  thresh.min <- 0
  pct.1 <- round(x = rowSums(x = object[features, cells.1, drop = FALSE] > thresh.min)/length(x = cells.1), digits = 3)
  pct.2 <- round(x = rowSums(x = object[features, cells.2,drop = FALSE] > thresh.min)/length(x = cells.2), digits = 3)
  data.1 <- mean.fxn(object[features, cells.1, drop = FALSE])
  data.2 <- mean.fxn(object[features, cells.2, drop = FALSE])
  fc <- (data.1 - data.2)
  fc.results <- as.data.frame(x = cbind(fc, pct.1, pct.2))
  colnames(fc.results) <- c(fc.name, "pct.1", "pct.2")
  return(fc.results)
}

#' Unknown function
PercentAbove <- function (x, threshold) 
{
  return(length(x = x[x > threshold])/length(x = x))
}

#' Unknown function
pos_pct <- function (x) 
{
  sum(x > 0)/length(x)
}
