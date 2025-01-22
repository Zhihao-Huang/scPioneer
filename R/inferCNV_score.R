#' Calculate CNV score
#' 
#' @export
cal_cnvScore <- function(infercnv_obj=NULL,genefile.path=NULL,
                         do.kmeans=FALSE,
                         k=7,do.correlation=FALSE,percent_top=5,
                         out.dir='.',out.prefix='cnvScore',
                         save.heatmap = F, 
                         min.cnv=0.4,max.cnv=1.2,
                         plot = F,
                         return.cordata = T){
  if (is.null(genefile.path)) {
    genefile.path <- "/storage2/hlinglab/jasper/project/tumor_evolution/result/100.infercnv/mm10_gene_file.txt"
  }
  geneFile <- read.table(genefile.path,header = F,sep = "\t",stringsAsFactors = F)
  geneFile <- geneFile[!grepl('chr24|chr25',geneFile$V2),]
  
  if (!dir.exists(out.dir)) {dir.create(out.dir)}
  out.dir <- sub('/$','',out.dir)
  
  obj_infercnv <- readRDS(infercnv_obj)
  ### extract observation index && filter referecnce cells
  cell.index.df <- c() # cell annotation 
  for(n in 1:length(obj_infercnv@observation_grouped_cell_indices)){
    temp <- obj_infercnv@observation_grouped_cell_indices[[n]]
    barcode <- colnames(obj_infercnv@expr.data[,temp,drop=F])
    temp.df <- data.frame(cell.idx=obj_infercnv@observation_grouped_cell_indices[[n]],
                          group=rep(names(obj_infercnv@observation_grouped_cell_indices[n]),length(obj_infercnv@observation_grouped_cell_indices[[n]])),
                          CB=barcode)
    cell.index.df <- rbind(cell.index.df,temp.df)
  }
  
  ############## calculate cnv score
  expr <- obj_infercnv@expr.data[,cell.index.df$cell.idx, drop=FALSE]
  expr2 <- expr-1
  expr2 <- expr2 ^ 2
  CNV_score <- as.data.frame(colMeans(expr2))
  colnames(CNV_score) <- "CNV_score"
  CNV_score$CB <- rownames(CNV_score)
  if(min.cnv > max(CNV_score$CNV_score)) stop('min.cnv is out of range.')
  if(max.cnv < min(CNV_score$CNV_score)) stop('max.cnv is out of range.')
  if(!do.kmeans){
    write.table(CNV_score, file = paste0(out.dir,'/',out.prefix,'.cnvScore.txt'), quote = FALSE, sep = '\t', row.names = T,col.names = T)
  }else{
    
    ############## do clustering (kmeans)
    expr3 <- expr2
    gn <- rownames(expr3)
    rownames(geneFile) <- geneFile$V1
    sub_geneFile <-  geneFile[intersect(gn,geneFile$V1),]
    expr3 <- expr3[intersect(gn,geneFile$V1),as.character(cell.index.df$CB)]
    
    set.seed(20210418)
    # k <- 7
    kmeans.result <- kmeans(t(expr3), k)
    kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
    kmeans_df$CB <- rownames(kmeans_df)
    kmeans_df <- kmeans_df%>%inner_join(cell.index.df,by="CB") # merge with annotation dataframe
    kmeans_df_s <- arrange(kmeans_df,kmeans_class) # sort 
    rownames(kmeans_df_s)=kmeans_df_s$CB
    # filter cluster contained only one cell
    pos <- table(kmeans_df_s$kmeans_class) == 1
    if (any(pos)) {
      filter_c <- names(table(kmeans_df_s$kmeans_class))[pos]
      message(paste0('Cluster ',unique(filter_c),' only contain 1 cell, ', rownames(kmeans_df_s)[kmeans_df_s$kmeans_class == filter_c]),', excluded.')
      kmeans_df_s <- kmeans_df_s[!kmeans_df_s$kmeans_class %in% filter_c, ]
    }
    # kmeans_df_s$CB=NULL
    kmeans_df_s$kmeans_class=as.factor(kmeans_df_s$kmeans_class) # transform kmeans clustersto factor
    # head(kmeans_df_s)
    kmeans_cnv <- dplyr::left_join(x = kmeans_df_s, y = CNV_score, by = 'CB')
    # output kmeans_cnv socre table 
    write.table(kmeans_cnv, file = paste0(out.dir,'/',out.prefix,'.cnvScore_kmeans_K',k,'.txt'), quote = FALSE, sep = '\t', row.names = T,col.names = T)
    
    if(do.correlation){
      ###################################### CNV correlation ###################################
      # percent_top <- 5 # top 5 cnv-score cells
      
      top_MS_cells <- arrange(kmeans_cnv, desc(CNV_score))[1:round(nrow(kmeans_cnv)*percent_top/100),]$CB  # Top 5%
      avg_cnvProfile <- rowMeans(expr2[,top_MS_cells,drop=F])
      
      ## calculate correlation : corr using 1 cell vs. top_MS_cells
      kmeans_cnv_cor <- kmeans_cnv
      for(i in 1:nrow(kmeans_cnv)){
        kmeans_cnv_cor$COR[i] <-  cor(expr2[,i, drop=FALSE], avg_cnvProfile)
      }
      write.table(kmeans_cnv_cor, file = paste0(out.dir,'/',out.prefix,'.cnvScore_kmeans_K',k,'_correlation.txt'), quote = FALSE, sep = '\t', row.names = T,col.names = T)
    }
    if (plot) {
      p_vio <- ggplot(kmeans_cnv,aes(x=kmeans_class,y=CNV_score,fill=kmeans_class)) + geom_violin() + scale_fill_d3()
      
      # add heatmap anno
      top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), 
                                                     labels = 1:22,labels_gp = gpar(cex = 1.5)))
      # color_v=RColorBrewer::brewer.pal(8, "Dark2")[1:k]  # num of clusters
      color_v <- pal_d3("category20")(20)[1:k] # color for clusters
      names(color_v) <- as.character(1:k)
      class_col <- rainbow(length(unique(cell.index.df$group))) # color for cell type
      names(class_col) <- unique(cell.index.df$group)
      left_anno <- rowAnnotation(df = kmeans_df_s[,c('kmeans_class','group')],
                                 col=list(group=class_col,kmeans_class=color_v))
      
      heatmap_mt <- t(expr)[rownames(kmeans_df_s),sub_geneFile$V1]
      ht <- Heatmap(heatmap_mt, #
                    col = colorRamp2(c(min.cnv,1,max.cnv), c("#377EB8","#F0F0F0","#E41A1C")), # heatmap color limits
                    cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
                    column_split = factor(sub_geneFile$V2, paste("chr",1:22,sep = "")), 
                    column_gap = unit(2, "mm"),
                    
                    heatmap_legend_param = list(title = "Modified expression",
                                                direction = "vertical",title_position = "leftcenter-rot",
                                                at=c(min.cnv,1,max.cnv),legend_height = unit(3, "cm")),
                    raster_device = "CairoPNG",
                    top_annotation = top_anno,left_annotation = left_anno, # add annotation
                    row_title = NULL,column_title = NULL)
      
      
      
      # save replot heatmap
      if (save.heatmap) {
        pdf(paste0(out.dir,'/',out.prefix,'.cnv.heatmap_K',k,'.pdf'),width = 15,height = 10)
        draw(ht, heatmap_legend_side = "right")
        dev.off()
      }
      # save cnv score violin plot
      pdf(paste0(out.dir,'/',out.prefix,'.cnvScore_kmeans_K',k,'.pdf'),width = 15,height = 10)
      print(p_vio)
      dev.off()
      # plot cor
      p_point <- ggplot(kmeans_cnv_cor,aes(x=CNV_score,y=COR,color=kmeans_class)) + geom_point() + scale_color_d3()
      # save cnv score-cnv correlation plot
      pdf(paste0(out.dir,'/',out.prefix,'.cnvScore_correlation.pdf'),width = 15,height = 10)
      print(p_point)
      dev.off()
    }
    if (return.cordata) {
      return(kmeans_cnv_cor)
    }
  }
}
