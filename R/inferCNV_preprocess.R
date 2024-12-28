#' Function of pre-processing data for inferCNV
#'
#' @param annotation Column of meta.data for distinguishing malignant and reference.
#' @param malignant_prefix Prefix of malignant names.
#' @param group Column of meta.data for annotation bar on the left side of CNV heatmap.
#' @examples 
#' datalist <- inferCNV_preprocess(object, meta.data,annotation, malignant_prefix = '^EP',
#' group,min.cell.per.group = 30,gene.anno.gtf,outdir,subset.malig = F,
#' subset.ref = T,subset.ref.prop = 0.5,seed = 123)
#' infercnv_obj = CreateInfercnvObject(raw_counts_matrix=datalist$matrix,
#' annotations_file=datalist$annotations_file,delim="\t",
#' gene_order_file=datalist$gene_order_file,
#' ref_group_names = datalist$ref_group_names)
#' infercnv_obj = infercnv::run(infercnv_obj,cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
#' out_dir=outdir, cluster_by_groups=TRUE, denoise=TRUE,HMM=TRUE)
#' @export
inferCNV_preprocess <- function(object, meta.data,
                                annotation, malignant_prefix = '^EP',
                                group,
                                min.cell.per.group = 30,
                                gene.anno.gtf,
                                outdir,
                                subset.malig = F,subset.malig.num = 10000,
                                subset.ref = T,subset.ref.prop = 0.5,
                                seed = 123) {
  EPIpos <- grepl(malignant_prefix, meta.data[,annotation])
  meta.data[,annotation][!EPIpos] <- meta.data[,annotation][!EPIpos]
  #EPI
  EPItypes <- unique(meta.data[,annotation][EPIpos])
  EPImeta <- meta.data[meta.data[,annotation] %in% EPItypes,]
  subcell_EPI <- EPImeta[,group, drop = F]
  subcell_EPI$Cellnames <- rownames(subcell_EPI)
  colnames(subcell_EPI) <- c('Cluster', 'Cellnames')
  if (subset.malig) {
    set.seed(seed)
    subcell_EPI <- as.data.frame(angrycell::subsetID(subcell_EPI,num = subset.malig.num))
    colnames(subcell_EPI)[1] <- 'Cluster'
  }
  # filter clusters whose number is less than 30.
  poscl <- table(subcell_EPI$Cluster) < min.cell.per.group
  selected_cl <- names(poscl)[!poscl]
  subcell_EPI <- subcell_EPI[subcell_EPI$Cluster %in% selected_cl,]
  dim(subcell_EPI)
  # ref cells
  identmat_NC <- meta.data[!EPIpos, group, drop = F]
  identmat_NC[,group] <- as.vector(identmat_NC[,group])
  identmat_NC$Cellnames <- rownames(identmat_NC)
  colnames(identmat_NC) <- c('Cluster', 'Cellnames')
  if (subset.ref) {
    set.seed(seed)
    num <-round(nrow(subcell_EPI) * subset.ref.prop,0)
    subcell_NC <- as.data.frame(angrycell::subsetID(identmat_NC,num = num))
    colnames(subcell_NC)[1] <- 'Cluster'
  }
  subcell <- rbind(subcell_NC, subcell_EPI)
  # subset object.
  subs <- subset(object, cells = subcell$Cellnames)
  #######################################################cell annotation
  subcell <- data.frame(subcell, stringsAsFactors = F)
  rownames(subcell) <- subcell$Cellnames
  anno <- cbind(rownames(subs@meta.data),as.vector(subcell[rownames(subs@meta.data),'Cluster']))
  colnames(anno) <- c('Cellname', 'Cluster')
  write.table(anno, file = paste0(outdir,'inferCNV_anno.txt'), col.names = F, row.names = F, quote=F, sep = '\t')
  #######################################################gene position
  gtf <- rtracklayer::import(gene.anno.gtf)
  gene_positions <- as.data.frame(gtf)
  # only gene.
  #gene_positions <- gene_positions[gene_positions$type == 'gene',]
  # remove duplicate
  gene_positions = gene_positions[!duplicated(gene_positions[, 'gene_name']), ]
  # change chromosome names.
  gene_positions$seqnames <- as.vector(gene_positions$seqnames)
  gene_positions[which(gene_positions[, 'seqnames'] == "X"), 'seqnames'] = 23
  gene_positions[which(gene_positions[, 'seqnames'] == "Y"), 'seqnames'] = 24
  gene_positions[which(gene_positions[, 'seqnames'] == "MT"), 'seqnames'] = 25
  ##other scaffolds = 26
  gene_positions[which(nchar(gene_positions[, 'seqnames']) > 2), 'seqnames'] = 26
  gene_pos = gene_positions[order(as.numeric(gene_positions[,'seqnames']), decreasing = F), ]
  head(gene_pos)
  gene_pos$seqnames <- paste0('chr',gene_pos$seqnames)
  gene_order_file <- gene_pos[,c('gene_name','seqnames','start','end')]
  colnames(gene_order_file) <- c('hgnc_symbol', 'chromosome_name', 'start_position',
                                 'end_position')
  ##common genes
  commgene <- intersect(rownames(object), gene_order_file[,'hgnc_symbol'])
  # save gene position
  gene_order_file <- gene_order_file[gene_order_file[,'hgnc_symbol'] %in% commgene, ]
  write.table(gene_order_file, file = paste0(outdir,'inferCNV_gene_order_file.txt'), col.names = F,row.names = F, quote=F, sep = '\t')
  ############################################################raw counts
  matrix <- as.matrix(subs@assays$RNA@counts)
  matrix <- matrix[commgene,]
  ### return data list
  datalist <- list(raw_counts_matrix = matrix,
                   annotations_file = paste0(outdir,'inferCNV_anno.txt'),
                   gene_order_file = paste0(outdir,'inferCNV_gene_order_file.txt'),
                   ref_group_names = unique(subcell_NC$Cluster),
                   annotations = anno,
                   gene_order = gene_order_file
  )
  return(datalist)
}
