#' List of all category in MsigDB.
#' 
#' @export
list_all_GSEA_categories <- function () {
    df <- data.frame(category = c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"), 
        sub_category = c("-", "-", "CGP,CP,BIOCARTA,KEGG,PID,REACTOME,WikiPathways", 
        "MIR,mkRDB,MIR_Legacy,TFT,GTRD,TFT_Legacy", "CGN,CM", 
        "GO,BP,CC,MF,HPO", "-", "ImmuneSigDB,VAX", "_"), genesets = c("hallmark gene sets", 
        "positional gene sets", "curated gene sets", "regulatory target gene sets", 
        "computational gene sets", "ontology gene sets", "oncogenic signature gene sets", 
        "immunologic signature gene sets", "cell type signature gene sets"), 
        stringsAsFactors = F)
    return(df)
}


#' Read gtm file from MSigDB and export data.frame for AngryCellDB.
#' @export
readGMT <- function(gmtfile){
  x <- readLines(gmtfile)
  res <- strsplit(x, "\t")
  names(res) <- sapply(res, function(x) x[1])
  res <- lapply(res, function(x) paste(x[-c(1:2)], collapse = ','))
  db <- stack(res)
  db <- db[, c("ind", "values")]
  colnames(db) <- c("Annotation", "geneSymbol")
  return(db)
}


#' Hypergeometirc test by phyper function
#' @param N number of set with all the genes to consider in the analysis.
#' @export
get_terms_P_MsigDB <- function(candidate, N, dbname) {
  pos <- dbname %in% names(MSigDBlist)
  if (sum(pos) < length(dbname)) {
    stop(paste0(paste(dbname[!pos],collapse = ','), ' not in MsigDB.'))
  }
  db <- MSigDBlist[dbname]
  db <- lapply(seq_along(db), function(x,y,i){cbind(x[[i]], y[[i]])},x = db, y = names(db))
  db <- do.call(rbind,db)
  colnames(db) <- c('Annotation','geneSymbol','Database')
  pvalue <- c()
  termpct <- c()
  Gene_in_Term <- c()
  for (t in 1:nrow(db)) {
    term_gene <- gsub(" ", "", strsplit(db$geneSymbol[t], ",")[[1]])
    m = length(term_gene)
    n = N - m
    k = length(candidate)
    x = sum(candidate %in% term_gene)
    termpct <- c(termpct, paste0(x, "/",m))
    Gene_in_Term_set <- intersect(candidate, term_gene)
    Gene_in_Term_set <- paste(Gene_in_Term_set, collapse = ",")
    Gene_in_Term <- c(Gene_in_Term, Gene_in_Term_set)
    p <-  phyper(q=x - 1, m=m, n=n, k=k, lower.tail=FALSE)
    if (p < 0) {p = 0}
    pvalue <- c(pvalue, p)
  }
  adjp <- p.adjust(pvalue, "BH")
  termp <- data.frame(as.vector(db$Database), as.vector(db$Annotation), pvalue, adjp, Gene_in_Term,
                      termpct, as.vector(db$geneSymbol), stringsAsFactors = F)
  termp <- termp[order(termp[, 3]), ]
  colnames(termp) <- c("Database","Description", "Term.p.value", "p.adjust", "Gene_in_Term",
                       "GeneRatio", "All_Term_Genes")
  return(termp)
}


#' Function for gene ID transform
#' @export
changeID <- function(genelist, species = c('Homosapiens','mm10')[1],
                     biomart = "ENSEMBL_MART_ENSEMBL",
                     host = "dec2016.archive.ensembl.org",
                     keyType = c('SYMBOL',"ENTREZID","ENSEMBL")[3],
                     toType = c('SYMBOL',"ENTREZID","ENSEMBL")[1]) {
  if (species == 'Homosapiens') {
    dataset <- 'hsapiens_gene_ensembl'
    SYMBOL <- 'hgnc_symbol'
  }else if (species == 'mm10') {
    dataset <- 'mmusculus_gene_ensembl'
    SYMBOL <- 'mgi_symbol'
  }else{
    stop('Only Homosapiens or mm10 available.')
  }
  keyTypelist <- list('SYMBOL' = SYMBOL,
                      'ENTREZID' = 'entrezgene_id',
                      'ENSEMBL' = 'ensembl_gene_id')
  toTypelist <- list('SYMBOL' = SYMBOL,
                     'ENTREZID' = 'entrezgene_id',
                     'ENSEMBL' = 'ensembl_gene_id')
  ensembl = biomaRt::useMart(biomart = biomart,
                             dataset = dataset, host = host)
  geneID <- biomaRt::getBM(attributes = c(keyTypelist[[keyType]], toTypelist[[toType]]),
                           filters = toTypelist[[toType]], values = genelist)
  geneID$pos <- which(genelist %in% geneID[,1])
  message(paste0('Change ID. ',nrow(geneID)*100 / length(genelist),'% genes left.'))
  return(geneID)
}

#' To get all the MSigDB names.
#' @export
listMSigDBnames <- function(){
  print(MSigDBnames)
}

#' Function for enrichment of MSigDB
#' @param N number of set with all the genes to consider in the analysis.
#' @param markerlist Vector. Candidate different expression gene set for enrichment.
#' @param DEGfile Dataframe. DEG file result from FindAllMarker function in Seurat.
#' @param species Select one species, Homosapiens or mm10.
#' @param keyType gene ID format.
#' @param use.top.marker Select top n marker genes for each cluster in DEG file.
#' @param filter.gene.pattern Filtering genes with some pattern, such as 'MT-'.
#' @param ont Select a database. Using listMSigDBnames() to list all databases' names.
#' @param pvalueCutoff P value cut-off.
#' @param qvalueCutoff Adjust P value. Only BH method in this version.
#' @export
enrichMSigDB <- function(
  markerlist = NULL,
  DEGfile = NULL,
  N,
  species = c('Homosapiens','mm10')[1],
  keyType = c('SYMBOL',"ENTREZID","ENSEMBL")[1],
  use.top.marker = 100,
  filter.gene.pattern = NULL,
  ont='c2.cgp.v7.1',
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2)
{
  if (!is.null(markerlist)) {
    if (species == 'Homosapiens' & keyType == 'SYMBOL') {
      markerlist <- markerlist
    }else{
      markermat <- changeID(markerlist, species = species, keyType = keyType, toType = 'SYMBOL')
      markerlist <- markermat[,2]
    }
    MSigmat <- get_terms_P_MsigDB(markerlist, N, ont)
    MSigmat <- MSigmat[MSigmat$Term.p.value < pvalueCutoff & MSigmat$p.adjust < qvalueCutoff,]
    return(MSigmat)
  }else if (!is.null(DEGfile)) {
    #change gene ID if necessary
    if (species == 'Homosapiens' & keyType == 'SYMBOL') {
      DEGfile <- DEGfile
    }else{
      markermat <- changeID(DEGfile$gene, species = species, keyType = keyType, toType = 'SYMBOL')
      DEGfile <- DEGfile[markermat$pos,]
      DEGfile$gene <- markermat[,2]
    }
    #filter.gene.pattern
    if (!is.null(filter.gene.pattern)) {
      topn <- DEGfile %>%  filter(!grepl(filter.gene.pattern,DEGfile$gene)) %>%
        group_by(cluster) %>%  top_n(n =use.top.marker, wt = avg_logFC)
    }else{
      topn <- DEGfile %>%  group_by(cluster) %>% top_n(n =use.top.marker, wt = avg_logFC)
    }
    #enrichment for each cluster
    MSiglist <- list()
    for (t in unique(topn$cluster)){
      markers <- topn$gene[topn$cluster == t]
      MSienrichmat <-  get_terms_P_MsigDB(markers, N, ont)
      df <- as.data.frame(MSienrichmat)
      df$celltype <- t
      MSiglist[[paste0('cluster_',t)]] <- df
    }
    MSigmat <- do.call(rbind,MSiglist)
    MSigmat <- MSigmat[MSigmat$Term.p.value < pvalueCutoff & MSigmat$p.adjust < qvalueCutoff,]
    return(MSigmat)
  }else{
    stop('Please input gene set.')
  }
}


#' GO enrichment.
#' GO enrichment by clusterProfiler::enrichGO.
#'
#' @param markerlist A DEG dataframe from Seurat::FindAllMarkers.
#' @param use.top.marker Number of top markers for GO enrichment. Default is 100.
#' @param filter.gene.pattern Filtering gene.
#' @param ont Ontology. Parameter in clusterProfiler::enrichGO.
#' @param pAdjustMethod Parameter in clusterProfiler::enrichGO.
#' @param pvalueCutoff Parameter in clusterProfiler::enrichGO.
#' @param qvalueCutoff Parameter in clusterProfiler::enrichGO.
#' @param keyType Parameter in clusterProfiler::bitr.
#' @param ... Parameters in enrichGO.
#' @examples
#' markerlist <- FindAllMarkers(pbmc,only.pos = TRUE,logfc.threshold = 0.5)
#' GOmat <- clusterGO(markerlist)
#' @export
clusterGO <- function(
  markerlist,
  use.top.marker = 100,
  filter.gene.pattern = NULL,
  ont='BP',
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  keyType = 'SYMBOL',
  OrgDb = "org.Hs.eg.db",
  n.cores = 1,
  ...){
  if (!is.null(filter.gene.pattern)) {
    topn <- markerlist %>%  filter(!grepl(filter.gene.pattern,markerlist$gene)) %>%
      group_by(cluster) %>%  top_n(n =use.top.marker, wt = avg_logFC)
  }else{
    topn <- markerlist %>%  group_by(cluster) %>% top_n(n =use.top.marker, wt = avg_logFC)
  }
  GOlist <- parallel::mclapply(unique(topn$cluster), function(t) {
    markers <- topn$gene[topn$cluster == t]
    eg <- clusterProfiler::bitr(markers, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"),
                                OrgDb=OrgDb)
    genelist <- eg[,keyType]
    cluster0_markers_go <- clusterProfiler::enrichGO(genelist, 
                                                     OrgDb = OrgDb,
                                                     ont=ont,
                                                     pAdjustMethod = pAdjustMethod,
                                                     pvalueCutoff = pvalueCutoff,
                                                     qvalueCutoff = qvalueCutoff, 
                                                     keyType = keyType)
    df <- as.data.frame(cluster0_markers_go)
    if (nrow(df) == 0) {
      message(paste0('Warning: NO significant GO term found for cluster: ',t))
      df <- as.data.frame(matrix(nrow = nrow(df), ncol = ncol(df) + 1, 
                                dimnames = list(NULL,c(colnames(df),'celltype')))
                          )
    }else{
      df$celltype <- t
    }
    df
  }, mc.cores = n.cores)
  names(GOlist) <- unique(topn$cluster)
  GOmat <- do.call(rbind,GOlist)
  return(GOmat)
}

#' KEGG enrichment. KEGG enrichment by clusterProfiler::enrichKEGG.
#' 
#' @param markerlist A DEG dataframe from Seurat::FindAllMarkers.
#' @param use.top.marker Number of top markers for KEGG enrichment. Default is 100.
#' @param filter.gene.pattern Filtering gene.
#' @param pAdjustMethod	 Parameter in clusterProfiler::enrichKEGG.
#' @param pvalueCutoff Parameter in clusterProfiler::enrichKEGG.
#' @param qvalueCutoff Parameter in clusterProfiler::enrichKEGG.
#' @param ... Parameters in enrichKEGG.
#' 
#' @examples 
#' markerlist <- FindAllMarkers(pbmc,only.pos = TRUE,logfc.threshold = 0.5)
#' KEGGmat <- clusterKEGG(markerlist)
#' 
#' @export
clusterKEGG <- function (markerlist, 
                         use.top.marker = 100, 
                         organism = "hsa",
                         keyType = "kegg",
                         filter.gene.pattern = NULL, 
                         pAdjustMethod = "BH", 
                         pvalueCutoff = 0.05, qvalueCutoff = 0.2, 
          ...) {
  if (!is.null(filter.gene.pattern)) {
    topn <- markerlist %>% filter(!grepl(filter.gene.pattern,markerlist$gene)) %>%
      group_by(cluster) %>% top_n(n = use.top.marker, wt = avg_logFC)
  }
  else {
    topn <- markerlist %>% group_by(cluster) %>% top_n(n = use.top.marker, 
                                                       wt = avg_logFC)
  }
  KEGGlist <- list()
  for (t in unique(topn$cluster)) {
    markers <- topn$gene[topn$cluster == t]
    if (organism == 'hsa') {
      eg <- clusterProfiler::bitr(markers, fromType = "SYMBOL", 
                                  toType = c("ENTREZID", "ENSEMBL"), OrgDb = "org.Hs.eg.db")
      genelist <- eg[, "ENTREZID"]
    }else if (organism == 'mmu') {
      eg <- clusterProfiler::bitr(markers, fromType = "SYMBOL", 
                                  toType = c("ENTREZID", "ENSEMBL"), OrgDb = "org.Mm.eg.db")
      genelist <- eg[, "ENTREZID"]
    }else{stop('Only hsa or mmu allowed.')}
    cluster0_markers_KEGG <- clusterProfiler::enrichKEGG(genelist, 
                                                         organism = organism,
                                                         keyType = keyType,
                                                         pAdjustMethod = pAdjustMethod, 
                                                         pvalueCutoff = pvalueCutoff, 
                                                         qvalueCutoff = qvalueCutoff)
    df <- as.data.frame(cluster0_markers_KEGG)
    if (nrow(df) == 0) {
      message(paste0('Warning: NO significant KEGG term found for cluster: ',t))
      df <- as.data.frame(matrix(nrow = nrow(df), ncol = ncol(df) + 1, 
                                 dimnames = list(NULL,c(colnames(df),'celltype')))
      )
    }else{
      df$celltype <- t
    }
    KEGGlist[[t]] <- df
  }
  KEGGmat <- do.call(rbind, KEGGlist)
  return(KEGGmat)
}

#' Pathway analysis of DEGlist from FindMarkers using fGSEA.
#' 
#' pathway – name of the pathway;
#' pval – an enrichment p-value; 
#' padj – a BH-adjusted p-value;
#' log2err – the expected error for the standard deviation of the P-value logarithm.
#' ES – enrichment score, same as in Broad GSEA implementation; 
#' NES – enrichment score normalized to mean enrichment of random samples of the same size; 
#' size – size of the pathway after removing genes not present in 'names(stats)'. 
#' leadingEdge – vector with indexes of leading edge genes that drive the enrichment.
#' 
#' @param species Species of input data.
#' @param category Category of MsigDB, such as C2. Use list_all_GSEA_categories() to see all category. See https://www.gsea-msigdb.org/gsea/msigdb/ for details.
#' @param subcategory Sub-class of the category, such as KEGG.
#' @param simplify_pathways Simplify the pathways names. Lower the alphabets and remove '_'.
#' @param eps Minimum of P-value. See ?fgseaMultilevel.
#' @param nPermSimple Number of permutations for estimation of P-values. See ?fgseaMultilevel.
#' @param DEGdf list of data frames from FindMarkers. DEGdf should include 3 columns: avg_logFC(the second celltype is negative), cluster(2 celltypes) and gene. Recommend logfc.threshold = 0.
#' 
#' @examples 
#' DEG_GSEA <- parallel::mclapply(levels(Idents(pbmc)), function(x) {
#' mmk <- FindMarkers(pbmc, ident.1 = x, test.use = 'wilcox',
#' min.pct = 0.25,logfc.threshold = 0)
#' mmk$gene <- rownames(mmk)
#' mmk$cluster <- x
#' mmk$cluster[mmk$avg_logFC < 0] <- paste0(x,'_Others')
#' mmk
#' }, mc.cores = 10)
#' names(DEG_GSEA) <- levels(Idents(pbmc))
#' gsealist <- GSEA_all(DEG_GSEA, category = 'C2',subcategory = NULL, n.cores = 10)
#' fgsea_df <- do.call(rbind, gsealist)
#' 
#' @export
GSEA_all <- function (DEGdflist, species = "Homo sapiens", category = NULL, 
    subcategory = NULL, simplify_pathways = T, return.enrichment.plot = F, 
    enrichment.plot.pathways = NULL, eps = 0, nPermSimple = 1000, 
    n.cores = 4) 
{
    if (!is.null(subcategory)) {
        mdblist <- lapply(subcategory, function(x) {
            mdb <- msigdbr::msigdbr(species = species, category = category, 
                subcategory = x)
        })
        mdb <- do.call(rbind, mdblist)
    }
    else {
        mdb <- msigdbr::msigdbr(species = species, category = category)
    }
    fgsea_sets <- mdb %>% split(x = .$gene_symbol, f = .$gs_name)
    firstup <- function(x) {
        substr(x, 1, 1) <- toupper(substr(x, 1, 1))
        x
    }
    if (simplify_pathways) {
        paths <- gsub("^.*?_", "", names(fgsea_sets))
        paths <- tolower(paths)
        paths <- gsub("_", " ", paths)
        paths <- firstup(paths)
        names(fgsea_sets) <- paths
    }
    enrichplotlist <- list()
    clustername <- names(DEGdflist)
    gsealist <- pbapply::pblapply(clustername, cl = n.cores, 
        function(x) {
            DEGdf_one <- DEGdflist[[x]]
            cluster01.genes <- DEGdf_one %>% arrange(desc(avg_logFC)) %>% 
                dplyr::select(gene, avg_logFC)
            ranks <- tibble::deframe(cluster01.genes)
            if (return.enrichment.plot) {
                is.exist <- enrichment.plot.pathways %in% names(fgsea_sets)
                if (!is.exist) {
                  stop(paste0("Pathway '", enrichment.plot.pathways, 
                    "' was not in Category: ", category, " ", 
                    category_prefix))
                }
                p <- plotEnrichment(fgsea_sets[[enrichment.plot.pathways]], 
                  ranks) + labs(title = enrichment.plot.pathways)
                enrichplotlist[[x]] <- p
            }
            fgseaRes <- fgsea::fgseaMultilevel(pathways = fgsea_sets, 
                stats = ranks, eps = eps, nPermSimple = nPermSimple)
            fgseaRes <- fgseaRes[fgseaRes$NES > 0, ]
            fgseaRes$celltype <- x
            fgseaRes
        })
    names(gsealist) <- clustername
    if (return.enrichment.plot) {
        return(enrichplotlist)
    }
    return(gsealist)
}

#' Pathway analysis of DEG from FindMarkers using fGSEA. Only consider TWO cell types.
#' 
#' pathway – name of the pathway; pval – an enrichment p-value; padj – a BH-adjusted p-value; log2err – the expected error for the standard deviation of the P-value logarithm. ES – enrichment score, same as in Broad GSEA implementation; NES – enrichment score normalized to mean enrichment of random samples of the same size; size – size of the pathway after removing genes not present in 'names(stats)'. leadingEdge – vector with indexes of leading edge genes that drive the enrichment.
#' 
#' @param species Species of input data.
#' @param category Category of MsigDB, such as C2. Use list_all_GSEA_categories() to see all category. See https://www.gsea-msigdb.org/gsea/msigdb/ for details.
#' @param subcategory Sub-class of the category, such as KEGG.
#' @param simplify_pathways Simplify the pathways names. Lower the alphabets and remove '_'.
#' @param eps Minimum of P-value. See ?fgseaMultilevel.
#' @param nPermSimple Number of permutations for estimation of P-values. See ?fgseaMultilevel.
#' @param DEGdf list of data frames from FindMarkers. DEGdf should include 3 columns: avg_logFC(the second celltype is negative), cluster(2 celltypes) and gene. Recommend logfc.threshold = 0.
#' 
#' @examples 
#' celltype_1 <- 'CD8 T'
#' celltype_2 <- 'NK'
#' DEG_pair <- FindMarkers(pbmc, ident.1 = celltype_1, ident.2 = celltype_2, min.pct = 0.25, logfc.threshold = 0)
#' DEG_pair$gene <- rownames(DEG_pair)
#' DEG_pair$cluster <- celltype_1
#' DEG_pair$cluster[DEG_pair$avg_logFC < 0] <- celltype_2
#' GSEA_KEGG <- GSEA_pair(DEG_pair, category = 'C2', subcategory = 'KEGG')
#' GSEA_KEGG <- GSEA_KEGG[GSEA_KEGG$padj < 0.05, ]
#' GSEA_C2 <- GSEA_pair(DEG_pair, category = 'C2', subcategory = NULL)
#' GSEA_C2 <- GSEA_C2[GSEA_C2$padj < 0.05, ]
#' 
#' @export
GSEA_pair <- function (DEGdf, species = "Homo sapiens", category = NULL, subcategory = NULL, 
          simplify_pathways = T, return.enrichment.plot = F, enrichment.plot.pathways = NULL, 
          eps = 0, nPermSimple = 1000) 
{
  if (!is.null(subcategory)) {
    mdblist <- lapply(subcategory, function(x) {
      mdb <- msigdbr::msigdbr(species = species, category = category, 
                              subcategory = x)
    })
    mdb <- do.call(rbind, mdblist)
  }
  else {
    mdb <- msigdbr::msigdbr(species = species, category = category)
  }
  fgsea_sets <- mdb %>% split(x = .$gene_symbol, f = .$gs_name)
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  if (simplify_pathways) {
    paths <- gsub("^.*?_", "", names(fgsea_sets))
    paths <- tolower(paths)
    paths <- gsub("_", " ", paths)
    paths <- firstup(paths)
    names(fgsea_sets) <- paths
  }
  cluster01.genes <- DEGdf %>% arrange(desc(avg_logFC)) %>% 
    dplyr::select(gene, avg_logFC)
  ranks <- tibble::deframe(cluster01.genes)
  if (return.enrichment.plot) {
    is.exist <- enrichment.plot.pathways %in% names(fgsea_sets)
    if (!is.exist) {
      stop(paste0("Pathway '", enrichment.plot.pathways, 
                  "' was not in Category: ", category, " ", category_prefix))
    }
    p <- plotEnrichment(fgsea_sets[[enrichment.plot.pathways]], 
                        ranks) + labs(title = enrichment.plot.pathways)
    return(p)
  }
  cellname_1 <- unique(DEGdf$cluster[DEGdf$avg_logFC > 0])
  cellname_2 <- unique(DEGdf$cluster[DEGdf$avg_logFC < 0])
  message(paste0('Celltype 1: ',cellname_1))
  message(paste0('Celltype 2: ',cellname_2))
  rank1 <- ranks
  fgseaRes_cell1 <- fgsea::fgseaMultilevel(pathways = fgsea_sets, 
                                           stats = rank1, eps = eps, nPermSimple = nPermSimple)
  fgseaRes_cell1 <- fgseaRes_cell1[fgseaRes_cell1$NES > 0, 
  ]
  fgseaRes_cell1$celltype <- cellname_1
  rank2 <- -ranks
  fgseaRes_cell2 <- fgsea::fgseaMultilevel(pathways = fgsea_sets, 
                                           stats = rank2, eps = eps, nPermSimple = nPermSimple)
  fgseaRes_cell2 <- fgseaRes_cell2[fgseaRes_cell2$NES > 0, 
  ]
  fgseaRes_cell2$NES <- -fgseaRes_cell2$NES
  fgseaRes_cell2$celltype <- cellname_2
  fgseaRes <- rbind(fgseaRes_cell1, fgseaRes_cell2)
  # list to character
  fgseaRes$leadingEdge <- unlist(lapply(fgseaRes$leadingEdge, function(x) paste(x, collapse = ',')))
  return(fgseaRes)
}

