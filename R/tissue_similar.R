#' Preprocess Seurat object. Running PCA based on tissue specific HVGs.
#' 
#' @export
preprocessing_obj <- function(object, tissue_colnames, query_group, 
                              assay = 'RNA', npcs = 50, 
                              hvg.method = 'vst', hvg.nfeatures = 2000,
                              vars.to.regress = NULL) {
  # 2 groups: query tissue and other tissues
  group_pos <- object@meta.data[,tissue_colnames] == query_group
  other_pos <- object@meta.data[,tissue_colnames] != query_group
  message('Cell number of query tissue ',query_group,': ', sum(group_pos),
          '; Cell number of other tissues: ', sum(other_pos))
  message('Get HVGs from each tissue.')
  # get HVGs from each tissue.
  tissues <- as.vector(unique((object@meta.data[,tissue_colnames])))
  tissuelist <- lapply(tissues, function(x) {
    cellID <- colnames(object)[object@meta.data[,tissue_colnames] == x]
    subset(object, cells = cellID)
  })
  names(tissuelist) <- tissues
  HVGlist <- lapply(tissuelist, function(x) {
    x <- FindVariableFeatures(x, selection.method = hvg.method, 
                              nfeatures = hvg.nfeatures)
    x@assays$RNA@var.features
  })
  names(HVGlist) <- tissues
  # merge HVGs from each tissues.
  HVGs <- unique(unlist(HVGlist))
  message(paste0('Total number of HVGs: ',length(HVGs)))
  # scale data and run PCA
  object <- ScaleData(object = object, features = HVGs, assay = 'RNA',
                      vars.to.regress = vars.to.regress)
  object <- RunPCA(object = object, assay = 'RNA', features = HVGs, npcs = npcs)
  return(object)
}

#' Find knn among cells of all other tissues.
#' 
#' @export
run_queryknn <- function(pcamat, group_pos, other_pos, k = 30) {
  # split two groups in PCA space.
  mat_AS <- pcamat[group_pos,] 
  mat_other <- pcamat[other_pos,] 
  message('Running queryKNN...')
  # Running knn using default method KMKNN.
  time0 <- Sys.time()
  k.out <- BiocNeighbors::queryKNN(mat_other, mat_AS, k=k)
  time1 <- Sys.time()
  message(paste0('Costing ',round(time1-time0,0),' seconds.'))
  return(k.out)
}
#' knn network to matrix.
#' 
#' @export
get_nearest_mat <- function(k.out, tissue_labels, unique_labels) {
  # convert index to tissue labels.
  index <- k.out$index
  knn_labels <- sapply(1:ncol(index), function(x) tissue_labels[index[,x]])
  nearest_list <- lapply(1:nrow(index),  function(x) {
    freq_tab <- sort(table(knn_labels[x,]), decreasing = T)
    max_num <- max(freq_tab)
    most_freq_labels <- names(freq_tab)[freq_tab == max_num]
    # If there are more than 1 tissues shared common max number of knn, select one with shortest distance.
    if (length(most_freq_labels) > 1) {
      df <- data.frame(distance = k.out$distance[x,],
                       tissue = knn_labels[x,],
                       stringsAsFactors = F)
      df <- df %>% group_by(tissue) %>% summarise(dist_total = sum(distance))
      most_freq_labels <- df$tissue[df$dist_total %in% min(df$dist_total)]
      # if total distance is still same, return all.
      if (most_freq_labels > 1) {
        warning(paste0('Query index ',x,' exist multi origin of tissues: ',paste(most_freq_labels,collapse = ',')))
      }
    }
    # binary matrix with tissues as col names, 1 indicate nearest neighbor tissue. 
    mat <- matrix(0, ncol = length(unique_labels))
    colnames(mat) <- unique_labels
    mat[1, unique_labels %in% most_freq_labels] <- 1
    mat
  })
  nearest_mat <- do.call(rbind, nearest_list)
  return(nearest_mat)
}
#' Permutation test by shuffling tissue labels.
#' 
#' @export
permutation_test_tissue <- function(object, tissue_colnames, query_group, 
                                    pcamat, k.out,
                                    k = 30,
                                    seed = 123,
                                    n.perm = 1000,
                                    pval.cutoff = 0.05,
                                    n.core = 4) {
  # 2 groups: query tissue and other tissues
  group_pos <- object@meta.data[,tissue_colnames] == query_group
  other_pos <- object@meta.data[,tissue_colnames] != query_group
  # split two groups in PCA space.
  mat_AS <- pcamat[group_pos,] 
  mat_other <- pcamat[other_pos,] 
  # permutation test, shuffling labels 1000 times
  other_num <- nrow(mat_other)
  set.seed(seed)
  sim <- replicate(n.perm, sample(other_num,other_num))
  # sim is a matrix recording iteration of sample;
  # Number of sim's columns is shuffling times.
  tissue_labels <- as.vector(object@meta.data[other_pos, tissue_colnames])
  # running permutation test
  message(paste0('Running permutation test. Shuffling ',n.perm,' times.'))
  nearest_mat_list <- parallel::mclapply(1:n.perm, function(x) {
    shuffle_labels <- tissue_labels[sim[,1]]
    unique_labels <- unique(tissue_labels)
    nearest_mat <- get_nearest_mat(k.out, shuffle_labels, unique_labels)
    return(nearest_mat)
  }, mc.cores = n.core)
  # sum matrics
  freqmat <- Reduce("+", nearest_mat_list)
  pvalmat <- 1 - freqmat / n.perm
  signif_pos <- pvalmat < pval.cutoff
  freq_pos <- freqmat > 0
  net_pos <- signif_pos * freq_pos
  # matrix to df
  net_pos <- as.data.frame(net_pos)
  pvalmat <- as.data.frame(pvalmat)
  rownames(net_pos) <- colnames(object)[group_pos]
  rownames(pvalmat) <- colnames(object)[group_pos]
  ## Add tissues that were totally disconnected, inferred from levels of tissues in meta.data.
  glevels <- levels(object@meta.data[,tissue_colnames])
  if (!is.null(glevels)) {
    disconnect_g <- !glevels %in% c(colnames(net_pos), query_group)
    if (any(disconnect_g)) {
      add_0_g <- glevels[disconnect_g]
      net_pos[,add_0_g] <- 0
      pvalmat[,add_0_g] <- 0
    }
  }
  datalist <- list(net_pos = net_pos, pvalmat = pvalmat)
  return(datalist)
}

#' To quantify the cell similarities among different tissues.
#' 
#' Query for each cell's nearest neighbor among cells of all other tissues 
#' in the low-dimensional space. Inspired by Fig2A in 
#' https://doi.org/10.1016/j.cell.2019.10.003.
#' 
#' @param object Seurat object.
#' @param tissue_colnames The colnumn name of seurat meta data store tissue labels.
#' @param query_group Which label is query tissue.
#' @param assay Surat assays.
#' @param npcs Parameter in Seurat::RunPCA.
#' @param hvg.method Select method for Seurat::FindVariableFeatures.
#' @param hvg.nfeatures Number of HVG.
#' @param vars.to.regress parameter in Seurat::ScaleData
#' @param k Number of nearest neighbor of each query cell.
#' @param seed Random seed for permutation test.
#' @param n.perm times of permutation test.
#' @param pval.cutoff P value cut off.
#' @param n.core Number of processes are applied when running permutation test.
#' 
#' @return A binary matrix. Columns is tissues, rows is query cells. 
#' 1 indicate connected.
#' 
#' @examples 
#' main <- readRDS("/BGFS1/projectdata/project_xinhua_wangxipeng/result/rds_HGSOC/HGSOC_main_rm_cellcycle.rds")
#' net_mat <- get_neareast_tissue(object = main,tissue_colnames = 'Group_abb',
#' query_group = 'AS', vars.to.regress = 'percent.mt', n.core = 48)
#' 
#' @export
get_neareast_tissue <- function(object, 
                                tissue_colnames, query_group, 
                                use.exist.reduction = NULL,
                                assay = 'RNA', npcs = 50, 
                                hvg.method = 'vst', hvg.nfeatures = 2000,
                                vars.to.regress = NULL,
                                k = 30,
                                seed = 123,
                                n.perm = 1000,
                                pval.cutoff = 0.05,
                                n.core = 4) {
  if (!is.null(use.exist.reduction)) {
    pcamat <- object@reductions[[use.exist.reduction]]@cell.embeddings[,1:npcs]
  }else{
    # merge tissue specific HVGs and run PCA.
    object <- preprocessing_obj(object = object, 
                                tissue_colnames = tissue_colnames, 
                                query_group = query_group, 
                                assay = assay, npcs = npcs, 
                                hvg.method = hvg.method,
                                hvg.nfeatures = hvg.nfeatures,
                                vars.to.regress = vars.to.regress)
    pcamat <- object@reductions$pca@cell.embeddings[,1:npcs]
  }
  # Running knn using default method KMKNN.
  group_pos <- object@meta.data[,tissue_colnames] == query_group
  other_pos <- object@meta.data[,tissue_colnames] != query_group
  if (sum(group_pos) == 0) {
    stop(paste0('Could not find ',query_group, ' in object meta data.'))
  } 
  k.out <- run_queryknn(pcamat, group_pos, other_pos, k = k)
  # permutation test, shuffling labels 1000 times
  net_pos_list <- permutation_test_tissue(object, tissue_colnames, query_group, 
                                          pcamat, k.out,
                                          k = k,
                                          seed = seed,
                                          n.perm = n.perm,
                                          pval.cutoff = pval.cutoff,
                                          n.core = n.core)
  return(net_pos_list)
}


#' circle plot of network to show tissue origin.
#' 
#' @param margin.par Margin of plot. Order: bottom, left, top, right.
#' 
#' @export
draw_circle <- function(net_pos,
                        meta,
                        group,
                        query.group.name,
                        annotation,
                        group.order = NULL,
                        cell.order = NULL,
                        color.sector = NULL,
                        color.link = NULL,
                        rev_origin = F,
                        frac.cutoff = 0,
                        only.max.tissue = F,
                        alpha.link = 0.4,
                        labels.y = 1.5,
                        legend.title = '',
                        legend.title.position = "topleft",
                        legend.direction = c("horizontal", "vertical" )[1],
                        legend.adjust.x = unit(4, "mm"),
                        legend.adjust.y = unit(4, "mm"),
                        legend.adjust = c("left", "bottom"),
                        legend.nrow = NULL,
                        legend.lwd = 2,
                        legend.fontsize = 10,
                        margin.par = c(1, 1, 1, 1),
                        track.height = 0.1,
                        return.df = F) {
  net_pos$Cellname <- rownames(net_pos)
  netpost <- melt(net_pos, id.vars = 'Cellname', variable.name = group)
  netpost <- inner_join(netpost, meta[,c('Cellname', annotation)])
  colnames(netpost) <- c("Cellname",  "Group", "value","Annotation")
  df <- netpost %>% group_by(Annotation,Group) %>% summarise(count = sum(value == 1)) 
  dft <- data.frame(table(meta[,group]))
  colnames(dft) <- c('Group','total_count')
  df <- left_join(df, dft)
  df$target_frac <- df$count / df$total_count
  df$query_frac <- df$count / sum(meta[,group] == query.group.name)
  
  if (!is.null(cell.order)) {
    df$Annotation <- factor(df$Annotation, levels = cell.order)
    df <- df[order(df$Annotation),]
  }
  Annotation <- as.vector(unique(df$Annotation))
  type_num <- length(Annotation)
  if (is.null(color.link)) {
    color.link <- scPalette1(type_num)
  }
  colordf <- data.frame(Annotation = Annotation,
                        color = color.link[1:type_num],
                        stringsAsFactors = F)
  df <- left_join(df, colordf)
  df <- df %>% group_by(Annotation) %>% mutate(type_sum = sum(count)) %>% 
    mutate(type_frac = count / type_sum)
  # target tissue with max proportion
  df <- df %>% group_by(Annotation) %>% mutate(is_max = count == max(count))
  if(return.df) {
    return(df)
  }
  
  # Create circle
  plot.new()
  circle_size <- unit(1, "snpc")
  pushViewport(viewport(x = 0.5, y = 1, width = circle_size, height = circle_size,
                        just = c("center", "top")))
  par(omi = gridOMI(), new = TRUE)
  if (is.null(group.order)) {
    sectors = as.vector(unique(meta[,group]))
  }else{
    sectors <- group.order
  }
  
  par(mar = margin.par)
  circos.par(cell.padding = c(0, 0, 0, 0))
  circos.initialize(sectors, xlim = c(0, 1))
  sec_num <- length(sectors)
  if (is.null(color.sector)) {
    color.sector <- scPalette2(sec_num)
  }
  circos.track(ylim = c(0, 1), track.height = track.height,
               bg.col = color.sector[1:sec_num], bg.border = NA)
  df$Group <- as.vector(df$Group)
  # add labels
  adj_sectors <- c(sectors[sec_num], sectors[-sec_num])
  labeldf <- data.frame(x = 1:sec_num - 0.5,
                        y = labels.y,
                        labels = adj_sectors,
                        stringsAsFactors = F)
  for (i in 1:nrow(labeldf)) {
    circos.text(labeldf$x[i],labeldf$y[i],labeldf$labels[i])
  }
  
  # add links
  # se is a 2-elements vector tell which two tissues linked;
  # se_1_frac_begin is a matrix which indicate position in query tissue where lines begin.
  # se_1_frac is a 2-elements vector indicate the range of link in query tissue.
  # se_2_frac is a 2-elements vector indicate the range of link in target tissue.
  se_1_frac_begin = 0
  se_2_frac_begin = matrix(0, ncol = length(sectors) - 1)
  colnames(se_2_frac_begin) <- sectors[-which(sectors %in% query.group.name)]
  for(i in 1:nrow(df)) {
    tissue = df$Group[i]
    se = c(query.group.name,tissue)
    se_1_frac_end = se_1_frac_begin + df$query_frac[i]
    se_2_frac_end = se_2_frac_begin[1,tissue] + df$target_frac[i]
    se_1_frac = c(se_1_frac_begin, se_1_frac_end)
    se_2_frac = c(se_2_frac_begin[1,tissue], se_2_frac_end)
    if (rev_origin) {
      se_2_frac <- 1 - se_2_frac
    }
    # filter links
    if (!is.null(frac.cutoff)) {
      if (df$type_frac[i] < frac.cutoff) {
        transparency = 1
      }else{
        transparency = alpha.link 
      }
    }
    if (only.max.tissue) {
      if (df$is_max[i]) {
        transparency = alpha.link
      }else{
        transparency = 1
      }
    }
    # plot links
    circos.link(se[1], se_1_frac, se[2], se_2_frac, 
                col = add_transparency(df$color[i], transparency), border = NA)
    se_1_frac_begin = se_1_frac_end
    se_2_frac_begin[1,tissue] = se_2_frac_end
  }
  
  # add legend 
  upViewport()
  #Annotation <- as.vector(unique(df$Annotation))
  lgd_num <- length(Annotation)
  link_colors <- unique(df$color)
  if (is.null(legend.nrow)) {
    legend.nrow <- lgd_num + 1
  }
  lgd_lines = Legend(at = Annotation, type = "lines", 
                     legend_gp = gpar(col = link_colors, lwd = legend.lwd),
                     labels_gp = gpar(fontsize = legend.fontsize),
                     nrow = legend.nrow, 
                     title_position = legend.title.position, 
                     title = legend.title)
  lgd_list = packLegend(lgd_lines, direction = legend.direction)
  draw(lgd_list, 
       x = legend.adjust.x, 
       y = legend.adjust.y,
       just = legend.adjust)
  circos.clear()
}


#' Preprocess net data, return all connect cells.
#' 
#' @export
preprocess_netdata <- function(net_pos,
                               meta,
                               group,
                               annotation,
                               samples) {
  net_pos$Cellname <- rownames(net_pos)
  netpost <- melt(net_pos, id.vars = 'Cellname', variable.name = group)
  netpost <- inner_join(netpost, meta[,c('Cellname', annotation, samples)])
  colnames(netpost) <- c("Cellname",  "Group", "value","Annotation",'Sample')
  df <- netpost[netpost$value == 1, -1]
  return(df)
}

#' statistics of meta data. return frequency.
#' 
#' @export
stat_meta_freq <- function(df) {
  df <- df[,c("Group", "Annotation", 'Sample')]
  df <- df %>% group_by(Sample, Annotation,Group) %>% mutate(count_g = n())
  df <- df %>% group_by(Sample, Annotation) %>% mutate(count_anno = n()) 
  df$frac <- df$count_g/df$count_anno
  df <- unique(df)
  df <- df %>% group_by(Annotation, Group) %>% mutate(ave = mean(frac))
  return(df)
}

#' Barplot faceted in row. Plot data is from stat_meta_freq().
#' 
#' @export
barplot_facetrow <- function(df,
                             color_palette = NULL,
                             test.method = c("t.test", "wilcox.test"),
                             adjust.method = c("BH","holm", "hochberg", 
                                               "hommel", "bonferroni", "BY", 
                                               "fdr", "none"),
                             p.label = c("p.adj.signif", "p.adj", "p"),
                             add.signif = T,
                             return.stat = F) {
  test.method <- match.arg(arg = NULL, choices = test.method)
  adjust.method <- match.arg(arg = NULL, choices = adjust.method)
  p.label <- match.arg(arg = NULL, choices = p.label)
  stat.test <- dfgg %>% group_by(Annotation) 
  if (test.method == 't.test') {
    stat.test <- stat.test %>% 
      rstatix::t_test(frac ~ Group) %>%
      rstatix::adjust_pvalue(method = adjust.method) %>%
      rstatix::add_significance("p.adj") 
  }
  if (test.method == 'wilcox.test') {
    stat.test <- stat.test %>% 
      rstatix::wilcox_test(frac ~ Group) %>%
      rstatix::adjust_pvalue(method = adjust.method) %>%
      rstatix::add_significance("p.adj") 
  }
  if (return.stat) {
    return(stat.test)
  }
  stat.test <- stat.test[stat.test$p.adj.signif != 'ns',]
  stat.test <- stat.test  %>% rstatix::add_xy_position()
  plotdata <- unique(dfgg[,c("Group", "ave","Annotation")])
  p <- ggplot() + 
    geom_bar(data = dfgg, aes(Group, ave, color = Group),
             stat = 'identity', position = 'dodge') + theme_bw()+ 
    facet_grid(cols = vars(Annotation), scales = "free") + 
    ylab('Fraction of ascites cells') + xlab('') +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  if (is.null(color_palette)) {
    color_palette = scPalette2(cell_num+1)[-4]
  }
  p <- p + scale_color_manual(values = color_palette)
  # stat_pvalue_manual only scale color, not fill.
  if (add.signif) {
    p <- p + ggpubr::stat_pvalue_manual(stat.test, label = p.label)
  }
  p <- p + geom_bar(data = dfgg, aes(Group, ave, fill = Group),
                    stat = 'identity', position = 'dodge')
  p <- p + scale_fill_manual(values = color_palette)
  return(p)
}










