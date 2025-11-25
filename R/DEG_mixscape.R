#' DEG using mixscape
#' 
#' @examples
#' eccite <- LoadData(ds = "thp1.eccite")
#' DEG <- DEG_mixscape(eccite, groupname = 'gene', nt.class.name = 'NT',target_gene = 'STAT2' )
#' 
#' @export
DEG_mixscape <- function(object, groupname, target_gene, nt.cell.class = "NT", 
                         covariates = NULL, ndims = NULL, num.neighbors = 20,
                         assay = 'RNA', new.assay.name = "PRTB",
                         split.by = NULL,
                         nt.class.name = "NT", 
                         min.de.genes = 5, 
                         iter.num = 10, 
                         de.assay = "RNA",  test.use = 'wilcox',
                         min.pct = 0.1, 
                         logfc.threshold = 0.01) {
  ### RNA-based clustering is driven by confounding sources of variation
  # Prepare RNA assay for dimensionality reduction: 
  # Normalize data, find variable features and scale data.
  DefaultAssay(object = object) <- assay
  object <- NormalizeData(object = object) %>% FindVariableFeatures() %>% ScaleData()
  
  # Run Principle Component Analysis (PCA) to reduce the dimensionality of the data.
  object <- RunPCA(object = object)
  if (is.null(ndims)) {
    ndims <- findPC::findPC(object@reductions$pca@stdev, 
                            number = 50,
                            method = 'all',
                            aggregate = 'voting')
  }
  # Calculate perturbation signature (PRTB).
  object<- CalcPerturbSig(
    object = object, 
    assay = assay, 
    slot = "data", 
    gd.class = groupname, 
    nt.cell.class = nt.cell.class, 
    reduction = "pca", 
    ndims = ndims, 
    num.neighbors = num.neighbors, 
    split.by = split.by, 
    new.assay.name = new.assay.name)
  
  # Prepare PRTB assay for dimensionality reduction: 
  # Normalize data, find variable features and center data.
  DefaultAssay(object = object) <- new.assay.name
  
  # Use variable features from RNA assay.
  VariableFeatures(object = object) <- VariableFeatures(object = object[[assay]])
  object <- ScaleData(object = object, do.scale = F, do.center = T)
  
  # Run PCA to reduce the dimensionality of the data.
  object <- RunPCA(object = object, reduction.key = 'prtbpca', reduction.name = 'prtbpca')
  
  # Mixscape identifies cells with no detectable perturbation
  # Here, we are assuming each target gene class is a mixture of two Gaussian distributions 
  # one representing the knockout (KO) and the other the non-perturbed (NP) cells.
  object$crispr <- nt.class.name
  object$crispr[object@meta.data[,groupname] != nt.class.name] <- target_gene
  object <- RunMixscape(
    object = object, 
    assay = "PRTB", 
    slot = "scale.data", 
    labels = "crispr", 
    nt.class.name = nt.class.name, 
    min.de.genes = min.de.genes, 
    iter.num = iter.num, 
    de.assay = assay, 
    verbose = F,
    prtb.type = "KO")
  # Run DE analysis and visualize results of DEG among NT, NP and KO
  DefaultAssay(object) <- de.assay
  Idents(object = object) <- object$mixscape_class
  print(table(Idents(object = object)))
  DEG <- FindAllMarkers(object = object,
                     only.pos = T,
                     min.pct = min.pct, latent.vars = covariates,
                     logfc.threshold = logfc.threshold, 
                     test.use = test.use)
  return(DEG)
}
