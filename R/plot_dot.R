#' dotplot group
#'
#' @examples
#' ov_PT <- subset(ov_obj, Group_abb == 'PT')
#' ov_MT <- subset(ov_obj, Group_abb == 'MT')
#' datalist_PT <- RunCellphoneDB(object = ov_PT,annotation = 'celltype',outdir = './test/PT')
#' datalist_MT <- RunCellphoneDB(object = ov_MT,annotation = 'celltype',outdir = './test/MT')
#' #file_mean <- read.table(file='./test/PT/means.txt', header=T, sep='\t', check.names=F, as.is=T,stringsAsFactors = F)
#' #file_pval <- read.table(file='./test/PT/pvalues.txt', header=T, sep='\t', check.names=F, as.is=T,stringsAsFactors = F)
#' #datalist_PT <- list(mean = file_mean, pval = file_pval)
#' #file_mean <- read.table(file='./test/MT/means.txt', header=T, sep='\t', check.names=F, as.is=T,stringsAsFactors = F)
#' #file_pval <- read.table(file='./test/MT/pvalues.txt', header=T, sep='\t', check.names=F, as.is=T,stringsAsFactors = F)
#' #datalist_MT <- list(mean = file_mean, pval = file_pval)
#' df_PT <- interaction_data(datalist_PT$mean, datalist_PT$pval,cellchat.evidence = T)
#' df_MT <- interaction_data(datalist_MT$mean, datalist_MT$pval,cellchat.evidence = T)
#' df_PT$organ <- 'PT'
#' df_MT$organ <- 'MT'
#' df <- as.data.frame(rbind(df_PT, df_MT))
#' chemo_pathways <- 'CCL|CXCL|XCR|CX3C'
#' df_DC <- df[df$ligand == 'DC'& grepl(chemo_pathways, df$pathway_name),]
#' dotplot_g(df_DC, strip.angle = 90)
#' gridlist <-dotplot_g(df_DC, clip = 'off', panel.spacing = 0, x_position = 'bottom')
#' grid::grid.draw(gridlist)
#'
#' @export
dotplot_g <- function(plotdata, highlight = NULL, cell.order = NULL,gene.order = NULL,text.size=8,
                      text.x.size=8,text.y.size=8,legend.title.size=15,legend.text.size=12,
                      color.use= rev(brewer.pal(11,"RdYlBu")),facet.toward = 'col',
                      scale.min = NA, scale.max = NA, scale.size = c(0,1.5,3),
                      x_position = 'top', y_position = 'left',
                      legend.position = 'right', strip.text.size = 8,strip.angle = 45,
                      strip.vjust = 0, strip.hjust = 0,
                      axis.x.position = NULL, axis.y.position = NULL,
                      strip.text.switch = NULL, panel.spacing = 0.5, clip = 'on',
                      max.pair = 50){
  if (!is.null(highlight)) {
    plotdata <- left_join(plotdata,as.data.frame(highlight),by='variable')
  }
  if (!is.na(x = scale.min)) {
    plotdata[plotdata$Expression < scale.min, "Expression"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    plotdata[plotdata$Expression > scale.max, "Expression"] <- scale.max
  }
  if(!is.null(cell.order)){plotdata$variable <- factor(plotdata$variable,levels = cell.order[cell.order %in% unique(as.vector(plotdata$variable))])}
  if(!is.null(gene.order)){plotdata$interacting_pair <- factor(plotdata$interacting_pair,levels = gene.order[gene.order %in% unique(as.vector(plotdata$interacting_pair))])}
  if (!'pval_group' %in% colnames(plotdata)) {
    plotdata$pval_group <- 'No significant'
    plotdata$pval_group[plotdata$pvalue < 0.05] <- '*0.05'
    plotdata$pval_group[plotdata$pvalue < 0.01] <- '*0.01'
  }
  plotdata$pval_group <- factor(plotdata$pval_group, levels = c('No significant', '*0.05', '*0.01'))
  if (length(unique(plotdata$interacting_pair)) > max.pair) {
    stop('There are too many interacting_pair-interacting_pair pairs.')
  }
  mytheme <- theme_bw() + theme(plot.title = element_blank(),#element_text(family="ArialMT", size=25, hjust = 0.5),
                                panel.grid= element_blank(),
                                axis.title = element_blank(),
                                axis.ticks.length = unit(0.3, "cm"),
                                axis.text=element_text(family="ArialMT", size = text.size),
                                axis.text.x=element_text(angle=60, vjust=0, hjust=0,size= text.x.size),
                                legend.title=element_text(family="ArialMT", size=legend.title.size),
                                legend.text=element_text(family="ArialMT", size=legend.text.size),
                                legend.background=element_rect(fill="transparent"),
                                legend.position = legend.position,
                                plot.margin = unit(c(1,1,1,5), "cm"))
  #mycol <- rev(brewer.pal(11,"RdYlBu"))
  p <- ggplot(plotdata,aes(x=organ, y=variable))
  if (!is.null(x_position)) {
    p <- p + scale_x_discrete(position = x_position)
  }
  if (!is.null(y_position)) {
    p <- p + scale_y_discrete(position = y_position)
  }
  if (!is.null(highlight)) {
    p <- p +  geom_segment(aes(x = variable, xend = variable, y = sig_interacting_pair),
                           yend = ifelse(x_position == 'top', Inf, 0),
                           color = "red", size = 0.2, show.legend = F)
    p <- p +  geom_segment(aes(y = sig_interacting_pair, yend = sig_interacting_pair, x = variable),
                           xend = ifelse(y_position == 'right', Inf, 0),
                           color = "red", size = 0.2, show.legend = F)
  }
  p <- p +  geom_point(aes(x=organ, y=variable, colour = value,
                           size = pval_group))
  p <- p + scale_size_manual(values = scale.size)
  p <- p + scale_color_gradientn(name='Average\nExpression', colors=color.use)
  #p <- p + new_scale_colour()
  if (!is.null(highlight)) {
    p <- p + ggnewscale::new_scale("color")
    p <- p +  geom_point(aes(x=organ, y=variable, color = 'red',
                             size = 5),shape = 21,show.legend = F)
  }
  p <- p +  guides(size=guide_legend(title=""))+mytheme
  if (facet.toward == 'row') {
    p <- p + facet_grid(rows = vars(interacting_pair), scales = "free",switch = strip.text.switch) +
      theme(legend.position = legend.position, panel.grid.major = element_line(colour = NA),
            panel.grid.minor = element_line(colour = NA)) +
      theme(axis.text.x = element_text(angle = 90, size = text.x.size),
            axis.text.y = element_text(size = text.y.size))
  }else if (facet.toward == 'col') {
    # p <- p + coord_flip()
    p <- p + facet_grid(cols = vars(interacting_pair), scales = "free",switch = strip.text.switch) +
      theme(legend.position = legend.position, panel.grid.major = element_line(colour = NA),
            panel.grid.minor = element_line(colour = NA)) +
      theme(axis.text.x = element_text(angle = 270, size = text.x.size),
            axis.text.y = element_text(size = text.y.size))
  }
  p <- p + xlab("") + ylab("")
  p <- p + theme(legend.title = element_text(size = legend.title.size),
                 strip.text = element_text(size = strip.text.size,angle = strip.angle,
                                           vjust=strip.vjust, hjust=strip.hjust),
                 strip.background = element_blank(),
                 strip.placement = 'outside',
                 panel.spacing = unit(panel.spacing, "lines"))
  if(!is.null(axis.x.position)) {
    p <- p + scale_x_discrete(position = axis.x.position)
  }
  if(!is.null(axis.y.position)) {
    p <- p + scale_y_discrete(position = axis.y.position)
  }
  #To avoid strip text cut off, set clip = 'off'. From https://stackoverflow.com/questions/49740215/ggplot-facet-grid-label-cut-off
  if (clip == 'off') {
    pg <- ggplotGrob(p)
    for(i in which(grepl("strip-", pg$layout$name))){
      pg$grobs[[i]]$layout$clip <- "off"
    }
    #grid::grid.draw(pg)
    return(pg)
  }else{
    return(p)
  }
}

#' hightlight dotplot
#'
#' @export
dotplot_multicolor <- function(exp,ident, highlight, cell.order = NULL,gene.order = NULL,
                               circle.size.as = c('expression percent','expressed cell number'),
                               text.size=8,text.x.size=8,
                               text.y.size=8,legend.title.size=8,legend.text.size=8,
                               do.scale.circle.size = T,
                               legend.key.size = 1,
                               legend.box= "horizonal",
                               coord_flip=T,
                               color.use=c('grey','yellow','red'),
                               noSig.circle.color = 'white',
                               gene_text_color = NULL,
                               scale.min = NA, scale.max = NA, x_position = 'bottom', y_position = 'left',
                               text.x.angle = 90, x.vjust = 1, x.hjust = 1,
                               text.y.angle = 0, y.vjust = 0.5, y.hjust = 0.5,
                               legend.title.fill = 'Average\nExpression',
                               legend.title.percent = "Percent of\nExpression",
                               legend.position = 'right', line.size = 0.5,
                               circle.border.size = 1.2,
                               plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")){
  # 20220827 start1
  circle.size.as <- match.arg(arg = NULL, choices = circle.size.as)
  # 20220827 end1
  exp <- exp[names(ident),,drop = F]
  exp$celltypes <- as.vector(ident)
  dat <- melt(exp, id.vars="celltypes")
  colnames(dat) <- c('celltypes','genes','expression')
  ave <- tapply(dat$expression,list(dat$celltypes,dat$genes),mean)
  pct <- tapply(dat$expression,list(dat$celltypes,dat$genes),function(x) sum(x > 0) / length(x) * 100 )
  # 20220827 start2
  if (circle.size.as == 'expressed cell number'){
    message('Circle size as number of expressed cells')
    cell_number <- as.data.frame(table(ident), row.names = 1)
    cell_number <- cell_number[rownames(pct),,drop = F]
    pct <- pct * cell_number$Freq / 100
  }
  # 20220827 end2
  plotdata_a <- melt(ave)
  plotdata_p <- melt(pct)
  plotdata <- data.frame(plotdata_a,plotdata_p$value)
  colnames(plotdata) <- c('Celltype','Gene','Expression','Percent_expr')
  colnames(highlight) <- c('Celltype','Gene','color')
  highlight$significant <- 'significant'
  plotdata <- left_join(plotdata,highlight,by=c('Celltype','Gene'))
  plotdata$color[is.na(plotdata$color)] <- noSig.circle.color
  plotdata$significant[is.na(plotdata$significant)] <- 'Not_significant'
  plotdata$significant2 <- plotdata$significant
  pos <- plotdata$significant2 == 'significant'
  plotdata$significant2[pos] <- paste0(plotdata$Celltype[pos],'_',plotdata$significant2[pos])
  if (!is.na(x = scale.min)) {
    plotdata[plotdata$Expression < scale.min, "Expression"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    plotdata[plotdata$Expression > scale.max, "Expression"] <- scale.max
  }
  if(is.null(cell.order)){
    cell.order <- unique(plotdata$Celltype)
  }
  plotdata$Celltype <- factor(plotdata$Celltype,levels = cell.order)
  if(!is.null(gene.order)){
    # 20240624 start1
   # plotdata$Gene <- factor(plotdata$Gene, levels = gene.order)
    p <- plotdata %>%
      mutate(Gene = factor(Gene, levels= rev(gene.order))) %>%
      arrange(Celltype,Gene) %>%
      ggplot(aes(x=Celltype, y=Gene))
    # 20240624 end1
  }else{
    gene.order <- unique(plotdata$Gene)
    plotdata$Gene <- factor(plotdata$Gene, levels = gene.order)
    p <- plotdata %>%
      arrange(Celltype) %>%
      ggplot(aes(x=Celltype, y=Gene))
  }
  
  line.color = p$data$color
  p <- p + geom_segment( aes(x=Celltype, xend=Celltype, y=Gene,
                             yend=ifelse(x_position == 'top', Inf, 0)),
                         color=line.color,
                         size=ifelse(p$data$significant == 'significant', line.size, NA) )
  p <- p + geom_segment( aes(x = Celltype, y = Gene, yend=Gene,
                             xend=ifelse(y_position == 'right', Inf, 0)),
                         color = line.color,
                         size=ifelse(p$data$significant == 'significant', line.size, NA) )
  p <- p +geom_point(aes(color = Expression, size = ifelse(Percent_expr==0, NA, Percent_expr)))
  if (!is.null(x_position)) {
    p <- p + scale_x_discrete(position = x_position)
  }
  if (!is.null(y_position)) {
    p <- p + scale_y_discrete(position = y_position)
  }
  # 20220827 start3
  if (do.scale.circle.size) {
    if (circle.size.as == 'expressed cell number'){
      breaklegend <- c(floor(max(pct)*1/4), floor(max(pct)*2/4) , floor(max(pct)*3/4), max(pct))
      p <- p + scale_size_continuous(range = c(-1, 5), breaks = breaklegend)
    }else{
      p <- p + scale_size_continuous(range = c(-1, 5), breaks = seq(25, 100, 25))
    }
  }
  # 20220827 end3
  p <- p + scale_color_gradientn(name= legend.title.fill, colors=color.use)
  #p <- p + new_scale_colour()
  p <- p + ggnewscale::new_scale("color")
  p <- p +  geom_point(aes(size = ifelse(Percent_expr==0, NA, Percent_expr)),
                       color = line.color, shape = 1,show.legend = F, stroke = circle.border.size)
  # 20240624 start2
  plotdata <- p$data
  # 20240624 end2
  plotdata$color_text <- gsub('white','black',plotdata$color)
  plotdata$significant <- factor(plotdata$significant, levels = c('significant','Not_significant'))
  
  if (is.null(gene_text_color)) {
    colordf <-plotdata[order(plotdata$Gene,plotdata$significant),c('Gene','color_text')]
    gene_text_color <- colordf[!duplicated(colordf$Gene),]$color_text
  }
  plotdata$Gene <- factor(plotdata$Gene, levels = rev(levels(plotdata$Gene)))
  colordf <-plotdata[order(plotdata$Celltype,plotdata$significant,plotdata$Gene),
                     c('Celltype','color_text')]
  cell_tick_color <- colordf[!duplicated(colordf$Celltype),]$color_text
  mytheme <- theme_bw() + theme(plot.title = element_blank(),#element_text(family="ArialMT", size=25, hjust = 0.5),
                                panel.grid= element_blank(),
                                axis.title = element_blank(),
                                axis.ticks.length = unit(0.2, "cm"),
                                axis.ticks.x = element_line(colour = cell_tick_color, size = 0.5),
                                axis.ticks.y = element_line(colour = gene_text_color, size = 0.5),
                                axis.text=element_text(family="ArialMT", size = text.size),
                                legend.title=element_text(family="ArialMT", size=legend.title.size),
                                legend.text=element_text(family="ArialMT", size=legend.text.size),
                                legend.background=element_rect(fill="transparent"),
                                legend.position = legend.position,
                                legend.key.size = unit(legend.key.size,'line'),
                                legend.box = legend.box,
                                axis.text.x = element_text(angle = text.x.angle, vjust = x.vjust, hjust = x.hjust, size = text.x.size),
                                axis.text.y = element_text(angle = text.y.angle, vjust = y.vjust, hjust = y.hjust,
                                                           color = gene_text_color),
                                plot.margin = plot.margin)
  p <- p + labs(size=legend.title.percent) + mytheme
  if(coord_flip){
    p <- p + coord_flip()
  }
  return(p)
}

#' DotPlot of all cell pairs.
#'
#' @param mfilePath File path to file mean.txt from result of cellphoneDB.
#' @param pfilePath File path to file pvalues.txt form result of cellphoneDB.
#' @param outdir Path to save dotplot PDF.
#' @param Rlib R labrary path which was needed by dotplot.
#' @param dotplot_func R script of dotplot.
#'
#' @examples
#' result_dir <- '/BGFS1/home/huangzh/workspace/project/r_package_development/biospiper/biospiper_V1.1.0/biospiper/test/'
#' dotplot_Rscript(mfilePath = paste0(result_dir, 'means.txt'), pfilePath = paste0(result_dir, 'pvalues.txt'),outdir = result_dir)
#' #dotplot PDF were saved in outdir/cellphoneDB/dotplot/.
#'
#' @export
dotplot_Rscript <- function(mfilePath, pfilePath, outdir,
                            Rlib = '/BGFS1/home/huangzh/workspace/project/single_cell_pipeline_development/BiosPipe_V0.1.0/lib/Rlib/BiosPipe_V0.1.0/',
                            dotplot_func = '/BGFS1/home/huangzh/workspace/project/single_cell_pipeline_development/BiosPipe_V0.1.0/biospipe/modules/cellphoneDB/cellphoneDB_dotplot.R') {
  mfile <- read.table(paste0(outdir_t, '/means.txt'), header=T, sep='\t', check.names=F, as.is=T,stringsAsFactors = F)
  pfile <- read.table(paste0(outdir_t, '/pvalues.txt'), header=T, sep='\t', check.names=F, as.is=T,stringsAsFactors = F)
  #save all dotplot in outdir/cellphoneDB/dotplot/
  system(paste0('Rscript ', dotplot_func, ' ',mfilePath, ' ',pfilePath, ' ', outdir, ' ', Rlib))
}

#' Merging two dotplots by gene-gene interaction lines.
#'
#' @param nordata A normalize matrix of gene expression.
#' @param ident A named vector/factor. The identity of cells. Could be Idents(object).
#' @param interdf A dataframe including 4 columns: ligandcell,receptorcell,ligandgene and receptorgene.
#' @param ligandcell Column name of interdf.
#' @param receptorcell Column name of interdf.
#' @param ligandgene Column name of interdf.
#' @param receptorgene Column name of interdf.
#' @param cell.order.li Cell order of ligand cells.
#' @param cell.order.re Cell order of receptor cells.
#' @param gene.order.li Gene order of ligand genes.
#' @param gene.order.re Gene order of receptor genes.
#' @param na.as Replace NA as numeric if NA exists.
#' @param width The relative widths to arrange 3 ggplot objects: heatmap + links + heatmap.
#' @param do.scale Scale gene expression to range of 0-1.
#'
#' @return ggplot object.
#'
#' @examples
#' library(reshape2)
#' library(ggplot2)
#' library(dplyr)
#' ov_PT <- subset(biospiper::ov_obj, Group_abb == 'PT')
#' datalist_PT <- RunCellphoneDB(object = ov_PT,annotation = 'celltype',outdir = './test/PT')
#' #file_mean <- read.table(file='./test/PT/means.txt', header=T, sep='\t', check.names=F, as.is=T,stringsAsFactors = F)
#' #file_pval <- read.table(file='./test/PT/pvalues.txt', header=T, sep='\t', check.names=F, as.is=T,stringsAsFactors = F)
#' #datalist_PT <- list(mean = file_mean, pval = file_pval)
#' df_PT <- interaction_data(datalist_PT$mean, datalist_PT$pval,cellchat.evidence = T)
#' chemo_li <- c(rownames(ov_PT)[grepl('^CCL',rownames(ov_PT))], rownames(ov_PT)[grepl('^CXCL',rownames(ov_PT))],'XCL1','XCL2','CX3CL1')
#' df_PT <- df_PT[df_PT$ligandgene %in% chemo_li, ]
#' df_PT <- df_PT[, c('ligand','receptor','ligandgene','receptorgene','pvalue')]
#' heatmapcolor <- c('#4493c2','white','#f5b08f','#b93032','#b93032','#b93032','#b93032')
#' color.use <- grDevices::colorRampPalette(heatmapcolor)(10)
#' #No highlight
#' p <- inter_dotplot(ov_PT@assays$RNA@data, Seurat::Idents(ov_PT), df_PT,
#' highlight= 'no', width = c(6,1,5), legend.box.li = 'vertical', color.use = color.use)
#' #highlight by single color
#' p <- inter_dotplot(ov_PT@assays$RNA@data, Seurat::Idents(ov_PT), df_PT,
#' highlight= 'singlecolor', width = c(6,1,5), legend.box.li = 'vertical')
#' #highlight by multicolor
#' p <- inter_dotplot(ov_PT@assays$RNA@data, Seurat::Idents(ov_PT), df_PT,
#' highlight= 'multicolor', width = c(6,1,5), legend.box.li = 'vertical'))
#'
#' @export
inter_dotplot <- function (nordata, ident,
                           interdf,
                           complex_input_path = NULL,
                           add.complex.gene.names = F,
                           highlight = c('no','multicolor','singlecolor'),
                           multicolor.line = NULL,
                           singlecolor.line.li = 'red',
                           singlecolor.line.re = 'red',
                           noSig.circle.color = 'white',
                           color.use = c('grey','#111111','black'),
                           p_thresh = 0.01,
                           ligandcell = 'ligand', receptorcell = 'receptor',
                           ligandgene = 'ligandgene', receptorgene = 'receptorgene',
                           cell.order.li = NULL, cell.order.re = NULL,
                           gene.order.li = NULL, gene.order.re = NULL,
                           circle.size.as = c('expression percent','expressed cell number'),
                           na.as = 0.05,
                           scale.min = NA, scale.max = NA,
                           do.scale.circle.size = T,
                           coord_flip = F,
                           line.size = 0.5,
                           text.x.angle = 45, x.vjust = 0, x.hjust = 0,
                           text.y.angle = 0, y.vjust = 0.5, y.hjust = 0.5,
                           text.size = 8,
                           text.y.size = 8, text.x.size = 10, legend.title.size = 10,
                           x_position.li = 'top',
                           x_position.re = 'top',
                           legend.text.size = 8,
                           legend.key.size = 1,
                           legend.box.li = "horizontal",
                           legend.box.re = "vertical",
                           legend.title.li = 'Ligand\nExpression',
                           legend.title.re = 'Receptor\nExpression',
                           legend.position.li = 'bottom',
                           legend.position.re = 'right',
                           width = c(4,1,4),
                           plot.margin.li = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
                           plot.margin.re = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
                           ...) {
  highlight <- match.arg(arg = NULL, choices = highlight)
  circle.size.as <- match.arg(arg = NULL, choices = circle.size.as)
  # add expression of complex into matrix
  complexdata <- biospiper::add_complex_to_matrix(nordata,
                                       add.complex.gene.names = add.complex.gene.names,
                                       complex_input_path = complex_input_path)
  if (add.complex.gene.names) {
    genedf <- data.frame(ligandgene = gsub(' \\(.*$','',rownames(complexdata)),
                         comb.name = rownames(complexdata),
                         stringsAsFactors = F)
    interdf <- left_join(interdf, genedf, by = 'ligandgene')
    compos <- !is.na(interdf$comb.name)
    # not to replace as complex name if complex name is same as symbol name existing in matrix.
    pos <- compos & (!interdf$ligandgene %in% rownames(nordata))
    interdf$ligandgene[pos] <- interdf$comb.name[pos]
    interdf$comb.name <- NULL
    colnames(genedf)[1] <- 'receptorgene'
    interdf <- left_join(interdf, genedf, by = 'receptorgene')
    compos <- !is.na(interdf$comb.name)
    pos <- compos & (!interdf$receptorgene %in% rownames(nordata))
    interdf$receptorgene[pos] <- interdf$comb.name[pos]
    interdf$comb.name <- NULL
  }
  nordata <- rbind(nordata, complexdata)
  
  ligene <- unique(interdf[,ligandgene])
  regene <- unique(interdf[,receptorgene])
  licell <- unique(interdf[,ligandcell])
  recell <- unique(interdf[,receptorcell])
  ident_li <- ident[ident %in% licell]
  ident_re <- ident[ident %in% recell]
  datali <- as.data.frame(Matrix::t(nordata[ligene,names(ident_li),drop = F]))
  datare <- as.data.frame(Matrix::t(nordata[regene,names(ident_re),drop = F]))
  highlight_li <- unique(interdf[interdf$pvalue < p_thresh,c(ligandcell, ligandgene)])
  highlight_re <- unique(interdf[interdf$pvalue < p_thresh,c(receptorcell, receptorgene)])
  colnames(highlight_li) <- c('Celltype','Gene')
  colnames(highlight_re) <- c('Celltype','Gene')
  if (highlight == 'multicolor') {
    if (is.null(multicolor.line)) {
      #ligand
      multicolor.line <- scPioneer::scPalette2(length(unique(highlight_li$Gene)))
      mcol <- data.frame(Gene = unique(highlight_li$Gene), color = multicolor.line, stringsAsFactors = F)
      highlight_li <- left_join(highlight_li, mcol, by = 'Gene')
      #receptor
      multicolor.line <- scPioneer::scPalette2(length(unique(highlight_re$Gene)))
      mcol <- data.frame(Gene = unique(highlight_re$Gene), color = multicolor.line, stringsAsFactors = F)
      highlight_re <- left_join(highlight_re, mcol, by = 'Gene')
    }
  }else if (highlight == 'singlecolor') {
    highlight_li$color <- singlecolor.line.li
    highlight_re$color <- singlecolor.line.re
  }else{
    noSig.circle.color <- singlecolor.line.re <- singlecolor.line.li <- 'black'
    highlight_li$color <- singlecolor.line.li
    highlight_re$color <- singlecolor.line.re
    line.size = 0
  }
  #left heatmap
  p1 <- dotplot_multicolor(datali,ident_li,highlight=highlight_li,
                           circle.size.as = circle.size.as,
                           color.use = color.use,
                           noSig.circle.color = noSig.circle.color,
                           x_position = x_position.li, y_position = 'right',coord_flip = F,
                           legend.position = legend.position.li,
                           cell.order = cell.order.li, gene.order = gene.order.li,
                           line.size = line.size,
                           text.x.angle = text.x.angle, x.vjust = x.vjust, x.hjust = x.hjust,
                           text.y.angle = text.y.angle, y.vjust = y.vjust, y.hjust = y.hjust,
                           legend.title.size = legend.title.size,
                           legend.text.size = legend.text.size,
                           legend.key.size = legend.key.size,
                           legend.box= legend.box.li,
                           legend.title.fill = legend.title.li,
                           do.scale.circle.size = do.scale.circle.size,
                           plot.margin = plot.margin.li,
                           ...)
  #right heatmap
  p2 <- dotplot_multicolor(datare,ident_re,highlight=highlight_re,
                           circle.size.as = circle.size.as,
                           color.use = color.use,
                           noSig.circle.color = noSig.circle.color,
                           y_position = 'left',x_position = x_position.re, coord_flip = F,
                           legend.position = legend.position.re,
                           cell.order = cell.order.re, gene.order = gene.order.re,
                           line.size = line.size,
                           text.x.angle = text.x.angle, x.vjust = x.vjust, x.hjust = x.hjust,
                           text.y.angle = text.y.angle, y.vjust = y.vjust, y.hjust = y.hjust,
                           legend.title.size = legend.title.size,
                           legend.text.size = legend.text.size,
                           legend.key.size = legend.key.size,
                           legend.box = legend.box.re,
                           legend.title.fill = legend.title.re,
                           do.scale.circle.size = do.scale.circle.size,
                           plot.margin = plot.margin.re,
                           ...)
  ##Prepare data for interaction lines.
  dfp <- interdf[interdf$pvalue < p_thresh,]
  interdata <- unique(interdf[, c(ligandgene,receptorgene,'pvalue')])
  interdata$inter_pairs <- paste0(interdata[,ligandgene],'_', interdata[,receptorgene])
  ligene_num <- length(unique(interdf[,ligandgene]))
  regene_num <- length(unique(interdf[,receptorgene]))
  ##adjust spaces to the bottom and the top of links.
  li_margin <- 1/(ligene_num * 2)
  re_margin <- 1/(regene_num * 2)
  #ligand data
  linkdata1 <- data.frame(x = rep(0,ligene_num),
                          ligandgene = unique(interdf[,ligandgene]))
  if (is.null(gene.order.li)) {
    gene.order.li <- ligene
  }
  # debug20240624. correction of gene order. reverse ligand and receptor order by rev(gene.order.li) and rev(gene.order.re).
  linkdata1$ligandgene <- factor(linkdata1$ligandgene, levels = rev(gene.order.li))
  linkdata1 <- linkdata1[order(linkdata1$ligandgene),]
  linkdata1$y.pos = seq(li_margin,1-li_margin, length.out = ligene_num)
  linkdata1 <- left_join(linkdata1, interdata, by = 'ligandgene')
  linkdata1 <- linkdata1[-which(colnames(linkdata1) == 'receptorgene')]
  colnames(linkdata1)[2] <- 'genename'
  #receptor data
  linkdata2 <- data.frame(x = rep(1,regene_num),
                          receptorgene = unique(interdf[,receptorgene]))
  if (is.null(gene.order.re)) {
    gene.order.re <- regene
  }
  linkdata2$receptorgene <- factor(linkdata2$receptorgene, levels = rev(gene.order.re))
  linkdata2 <- linkdata2[order(linkdata2$receptorgene),]
  linkdata2$y.pos <- seq(re_margin,1-re_margin, length.out = regene_num)
  linkdata2 <- left_join(linkdata2, interdata, by = 'receptorgene')
  linkdata2 <- linkdata2[-which(colnames(linkdata2) == 'ligandgene')]
  colnames(linkdata2)[2] <- 'genename'
  #
  linkdata <- rbind(linkdata1, linkdata2)
  
  #linkdata$inter_pairs[linkdata$pvalue >= p_thresh] <- NA
  p3 <-ggplot(linkdata, aes(x=x, y=y.pos)) +
    #geom_point(size = 2)+
    geom_line(data = linkdata[linkdata$pvalue < p_thresh,],
              aes(group = inter_pairs), size=0.5)+
    #geom_text(aes(label=genename),hjust=0.5, vjust=0.5) +
    #scale_color_manual(values = c("#4575B4","#D73027")) +
    theme_void()+
    expand_limits(x = c(0,1), y = c(0,1)) +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    theme(legend.position = 'top',
          plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))
  #patchwork
  p <- p1 + p3 + p2 + patchwork::plot_layout(widths = width)
  return(p)
}

#' dotplot
#'
#' @examples
#' datalist <- RunCellphoneDB(object = ov_obj,annotation = 'celltype',outdir = './test')
#' plist <- interaction_dotplot(datalist$mean, datalist$pval, ligandcell = 'DC',receptorcell.prefix = 'CD')
#' @export
interaction_dotplot <- function(file_m, file_p, ligandcell=NULL, receptorcell=NULL,
                                receptorcell.prefix = NULL, only.significant = T,
                                min.p.value = 0.05,
                                order.cells = NULL,
                                sort.receptorcell = F,
                                order.genepair.auto = F, object = NULL, annotation = NULL,
                                order.genepair = NULL,
                                max.cell = Inf,
                                merge.plot = F,
                                return.data = F,
                                save.pdf = F,
                                outdir.prefix = NULL,
                                scale.pt.size = c(0,3,6),
                                size.axis.text = 20, size.axis.text.x = 13,
                                size.legend.title=20, size.legend.text = 20,...){
  file_mean <- file_m[!duplicated(file_m$interacting_pair),]
  file_pval <- file_p[!duplicated(file_p$interacting_pair),]
  rownames(file_mean) <- file_mean$interacting_pair
  rownames(file_pval) <- file_pval$interacting_pair
  allreceptorcells <- gsub('^.*\\|','',colnames(file_mean)[12:ncol(file_mean)])
  allligandcells <- gsub('\\|.*$','',colnames(file_mean)[12:ncol(file_mean)])
  plotligandcells <- unique(allligandcells)
  plotreceptorcells <- unique(allreceptorcells)
  ##plot select cell types if assigning.
  if (!is.null(ligandcell)) {
    posli <- allligandcells %in% ligandcell
    if (sum(posli) == 0) {
      return(paste0('no ligand found in ',ligandcell))
    }
    plotligandcells <- ligandcell
  }
  if (!is.null(receptorcell)) {
    posre <- allreceptorcells %in% receptorcell
    if (sum(posre) == 0) {
      return(paste0('no receptor of found in ',receptorcell))
    }
    plotreceptorcells <- receptorcell
  }
  if (!is.null(receptorcell.prefix)) {
    posre <- grepl(receptorcell.prefix, allreceptorcells)
    if (sum(posre) == 0) {
      return(paste0('no receptor of found in ',receptorcell))
    }
    plotreceptorcells <- unique(allreceptorcells[posre])
  }
  plist <- list()
  for (li in plotligandcells) {
    re <- plotreceptorcells
    message(paste0('Dotplot of secret cell type: ',li, ' to the others.'))
    posli <- allligandcells %in% li
    posre <- allreceptorcells %in% re
    if (sum(posli) == 0 | sum(posre) == 0) {
      print('No significant pairs.')
      next
    }
    ##select ligand cell and receptor cell
    col_index <- c(2,which(posre&posli) + 11)
    selectfile_mean <- file_mean[, col_index]
    # data processing
    dat <- reshape2::melt(selectfile_mean, id.vars="interacting_pair",measure.vars = colnames(selectfile_mean)[-1])
    dat <- dat[!is.na(dat$value),]
    selectfile_pval <- file_pval[,c(2,12:ncol(file_pval))]
    datp <- reshape2::melt(selectfile_pval, id.vars="interacting_pair",measure.vars = colnames(selectfile_pval)[-1])
    colnames(datp) <- c(colnames(datp)[-3],'p_value')
    ##merge mean and pval
    dat$variable <- as.vector(dat$variable)
    datp$variable <- as.vector(datp$variable)
    datp <- datp[datp$variable %in% dat$variable,]
    suppressMessages(dat <- dat %>% inner_join(datp))
    dat$interacting_pair <- factor(dat$interacting_pair, levels=unique(dat$interacting_pair))
    if (nrow(dat) == 0) {
      print('No significant pairs.')
      next
    }
    #dat[dat$p_value < 0.0001, 'pval_group'] <- '**** 0.0001'
    #dat[dat$p_value >= 0.0001 & dat$p_value < 0.001, 'pval_group'] <- '*** 0.001'
    smallp <- dat$p_value >= 0.000 & dat$p_value < 0.01
    if (sum(smallp) != 0) {
      dat[smallp, 'pval_group'] <- '** 0.01'
    }
    midp <- dat$p_value >= 0.01
    if (sum(midp) != 0) {
      dat[midp , 'pval_group'] <- '* 0.05'
    }
    nosig <- dat$p_value > 0.05
    if (sum(nosig) != 0) {
      dat[nosig , 'pval_group'] <- 'No significant'
    }
    #dat$pval_group <- factor(dat$pval_group, levels = c('**** 0.0001', '*** 0.001', '** 0.01', '* 0.05'))
    dat$pval_group <- factor(dat$pval_group, levels = c('No significant','* 0.05', '** 0.01'))
    colnames(dat) <- c('interacting_pair','variable', 'means', 'p_value','pval_group')
    if (only.significant) {
      sign_pairs <- unique(dat$interacting_pair[dat$p_value < min.p.value])
      dat$interacting_pair <- as.vector(dat$interacting_pair)
      dat <- dat[dat$interacting_pair %in% sign_pairs,]
    }
    if (is.null(order.cells)) {
      dat$variable <- as.vector(dat$variable)
      order.cellpair <- unique(dat$variable)
    }else{
      order.cellpair <- order.cells
    }
    dat$receptorcell <- gsub('^.*\\|','',dat$variable)
    allcells <- unique(c(gsub('^.*\\|','',order.cellpair),li))
    if (sort.receptorcell) {
      allcells <- sort(allcells)
      dat$receptorcell <- factor(dat$receptorcell, levels = allcells)
      dat <- dat[order(dat$receptorcell),]
      order.cellpair <- unique(as.vector(dat$variable))
    }
    dat$variable <- factor(dat$variable, levels = order.cellpair)
    if (order.genepair.auto) {
      dat$receptorgene <- gsub('^.*_','',dat$interacting_pair)
      allreceptorgenes <- unique(dat$receptorgene)
      message('Auto-order gene pairs...')
      if (is.null(annotation)) {
        object$annotation <- Idents(object)
        annotation <- 'annotation'
      }
      order_re <- order_receptorgene(object,annotation,allreceptorgenes,allcells,...)
      dat$receptorgene <- factor(dat$receptorgene, levels = order_re)
      dat <- dat[order(dat$receptorgene),]
      dat$interacting_pair <- as.vector(dat$interacting_pair)
      order_pairs <- unique(dat$interacting_pair)
      dat$interacting_pair <- factor(dat$interacting_pair, levels = order_pairs)
    }
    if (!is.null(order.genepair)) {
      dat$interacting_pair <- factor(dat$interacting_pair, levels = order.genepair)
    }
    # plot theme
    mytheme <- theme_bw() + theme(plot.title = element_blank(),#element_text(family="ArialMT", size=25, hjust = 0.5),
                                  panel.grid= element_blank(),
                                  axis.title = element_blank(),
                                  axis.ticks.length = unit(0.3, "cm"),
                                  axis.text=element_text(family="ArialMT", size = size.axis.text),
                                  axis.text.x=element_text(angle=90, vjust=0.5, hjust=1,size=size.axis.text.x),
                                  legend.title=element_text(family="ArialMT", size=size.legend.title),
                                  legend.text=element_text(family="ArialMT", size=size.legend.text),
                                  legend.background=element_rect(fill="transparent"),
                                  plot.margin = unit(c(1,1,1,5), "cm"))
    mycol <- brewer.pal(11,"RdYlBu")
    ##seperate data to avoid  cons memory exhausted.
    if (length(unique(dat$variable)) > max.cell) {
      dat1 <- dat[dat$variable %in% unique(dat$variable)[1:max.cell],]
      dat2 <- dat[dat$variable %in% unique(dat$variable)[max.cell:length(unique(dat$variable))],]
      p1 <- ggplot(dat1, aes(x=variable, y=interacting_pair, colour = means, size=pval_group)) + geom_point() +
        scale_color_gradientn(name='Interaction\nScore', colors=rev(mycol)) +
        scale_size_manual(values = scale.pt.size) +
        coord_flip()+
        mytheme
      p2 <- ggplot(dat2, aes(x=variable, y=interacting_pair, colour = means, size=pval_group)) + geom_point() +
        scale_color_gradientn(name='Interaction\nScore', colors=rev(mycol)) +
        scale_size_manual(values = scale.pt.size) +
        coord_flip()+
        mytheme
      plist[[li]] <- p1+p2
      if (save.pdf) {
        liname <- paste0(unique(gsub('/','_',li)),collapse = '_and_')
        rename <- paste0(unique(gsub('/','_',re)),collapse = '_and_')
        outfile1 <- paste0(outdir.prefix,'_dotplot_',liname,'_To_other_',max.cell,'_receptor_celltypes.pdf')
        width1 <- 60 * length(unique(dat1$interacting_pair))/350 +20
        height1 <- ceiling(max.cell/5) + 6
        pdf(outfile1,width1,height1)
        print(p1)
        print(p2)
        dev.off()
      }
    }else{
      p <- ggplot(dat, aes(x=variable, y=interacting_pair, colour = means, size=pval_group)) + geom_point() +
        scale_color_gradientn(name='Interaction\nScore', colors=rev(mycol)) +
        scale_size_manual(values = scale.pt.size) +
        coord_flip()+
        mytheme
      plist[[li]] <- p
      if (save.pdf) {
        liname <- paste0(unique(gsub('/','_',li)),collapse = '_and_')
        rename <- paste0(unique(gsub('/','_',re)),collapse = '_and_')
        repnum <- length(unique(re))
        outfile <- paste0(outdir.prefix,'_dotplot_',liname,'_To_other_',repnum,'_receptor_celltypes.pdf')
        width <- 60 * length(unique(dat$interacting_pair))/350 +20
        height <- ceiling(repnum/5) + 6
        pdf(outfile,width,height)
        p <- ggplot(dat, aes(x=variable, y=interacting_pair, colour = means, size=pval_group)) + geom_point() +
          scale_color_gradientn(name='Interaction\nScore', colors=rev(mycol)) +
          coord_flip()+
          mytheme
        print(p)
        dev.off()
      }
    }
  }
  #return data
  if (return.data) {
    return(dat)
  }
  #return plot
  if (merge.plot) {
    p <- plist[[1]]
    p$data <- do.call(rbind, lapply(plist, function(x) x$data))
    return(p)
  }else{
    return(plist)
  }
}

#' dotplot
#'
#' @examples
#' datalist <- RunCellphoneDB(object = ov_obj,annotation = 'celltype',outdir = './test')
#' plist <- interaction_dotplot(datalist$mean, datalist$pval, ligandcell = 'DC',receptorcell.prefix = 'CD')
#' @export
interaction_dotplot2 <- function(file_m, file_p, ligandcell=NULL, receptorcell=NULL,
                                 receptorcell.prefix = NULL, only.significant = T,
                                 min.p.value = 0.05,
                                 order.cells = NULL,
                                 sort.receptorcell = F,
                                 order.genepair.auto = F, object = NULL, annotation = NULL,
                                 order.genepair = NULL,
                                 max.cell = Inf,
                                 return.data = F,
                                 save.pdf = NULL,
                                 scale.pt.size = c(0,3,6),
                                 size.axis.text = 20, size.axis.text.x = 13,
                                 size.legend.title=20, size.legend.text = 20,...){
  file_mean <- file_m[!duplicated(file_m$interacting_pair),]
  file_pval <- file_p[!duplicated(file_p$interacting_pair),]
  rownames(file_mean) <- file_mean$interacting_pair
  rownames(file_pval) <- file_pval$interacting_pair
  allreceptorcells <- gsub('^.*\\|','',colnames(file_mean)[12:ncol(file_mean)])
  allligandcells <- gsub('\\|.*$','',colnames(file_mean)[12:ncol(file_mean)])
  plotligandcells <- unique(allligandcells)
  plotreceptorcells <- unique(allreceptorcells)
  ##plot select cell types if assigning.
  if (!is.null(ligandcell)) {
    posli <- allligandcells %in% ligandcell
    if (sum(posli) == 0) {
      return(paste0('no ligand found in ',ligandcell))
    }
    plotligandcells <- ligandcell
  }
  if (!is.null(receptorcell)) {
    posre <- allreceptorcells %in% receptorcell
    if (sum(posre) == 0) {
      return(paste0('no receptor of found in ',receptorcell))
    }
    plotreceptorcells <- receptorcell
  }
  if (!is.null(receptorcell.prefix)) {
    posre <- grepl(receptorcell.prefix, allreceptorcells)
    if (sum(posre) == 0) {
      return(paste0('no receptor of found in ',receptorcell))
    }
    plotreceptorcells <- unique(allreceptorcells[posre])
  }
  #plist <- list()
  #for (li in plotligandcells) {
  re <- plotreceptorcells
  li <- plotligandcells
  message(paste0('Ligand cell type: ',paste(li, collapse = ', ')))
  message('')
  message(paste0('Target cell type: ',paste(re, collapse = ', ')))
  posli <- allligandcells %in% li
  posre <- allreceptorcells %in% re
  if (sum(posli) == 0 | sum(posre) == 0) {
    print('No significant pairs.')
    next
  }
  ##select ligand cell and receptor cell
  col_index <- c(2,which(posre&posli) + 11)
  selectfile_mean <- file_mean[, col_index]
  # data processing
  dat <- reshape2::melt(selectfile_mean, id.vars="interacting_pair",measure.vars = colnames(selectfile_mean)[-1])
  dat <- dat[!is.na(dat$value),]
  selectfile_pval <- file_pval[,c(2,12:ncol(file_pval))]
  datp <- reshape2::melt(selectfile_pval, id.vars="interacting_pair",measure.vars = colnames(selectfile_pval)[-1])
  colnames(datp) <- c(colnames(datp)[-3],'p_value')
  ##merge mean and pval
  dat$variable <- as.vector(dat$variable)
  datp$variable <- as.vector(datp$variable)
  datp <- datp[datp$variable %in% dat$variable,]
  suppressMessages(dat <- dat %>% inner_join(datp))
  dat$interacting_pair <- factor(dat$interacting_pair, levels=unique(dat$interacting_pair))
  if (nrow(dat) == 0) {
    print('No significant pairs.')
    next
  }
  #dat[dat$p_value < 0.0001, 'pval_group'] <- '**** 0.0001'
  #dat[dat$p_value >= 0.0001 & dat$p_value < 0.001, 'pval_group'] <- '*** 0.001'
  smallp <- dat$p_value >= 0.000 & dat$p_value < 0.01
  if (sum(smallp) != 0) {
    dat[smallp, 'pval_group'] <- '** 0.01'
  }
  midp <- dat$p_value >= 0.01
  if (sum(midp) != 0) {
    dat[midp , 'pval_group'] <- '* 0.05'
  }
  nosig <- dat$p_value > 0.05
  if (sum(nosig) != 0) {
    dat[nosig , 'pval_group'] <- 'No significant'
  }
  #dat$pval_group <- factor(dat$pval_group, levels = c('**** 0.0001', '*** 0.001', '** 0.01', '* 0.05'))
  dat$pval_group <- factor(dat$pval_group, levels = c('No significant','* 0.05', '** 0.01'))
  colnames(dat) <- c('interacting_pair','variable', 'means', 'p_value','pval_group')
  if (only.significant) {
    sign_pairs <- unique(dat$interacting_pair[dat$p_value < min.p.value])
    dat$interacting_pair <- as.vector(dat$interacting_pair)
    dat <- dat[dat$interacting_pair %in% sign_pairs,]
  }
  if (is.null(order.cells)) {
    dat$variable <- as.vector(dat$variable)
    order.cellpair <- unique(dat$variable)
  }else{
    order.cellpair <- order.cells
  }
  dat$receptorcell <- gsub('^.*\\|','',dat$variable)
  allcells <- unique(c(gsub('^.*\\|','',order.cellpair),li))
  if (sort.receptorcell) {
    allcells <- sort(allcells)
    dat$receptorcell <- factor(dat$receptorcell, levels = allcells)
    dat <- dat[order(dat$receptorcell),]
    order.cellpair <- unique(as.vector(dat$variable))
  }
  dat$variable <- factor(dat$variable, levels = order.cellpair)
  if (order.genepair.auto) {
    dat$receptorgene <- gsub('^.*_','',dat$interacting_pair)
    allreceptorgenes <- unique(dat$receptorgene)
    message('Auto-order gene pairs...')
    if (is.null(annotation)) {
      object$annotation <- Idents(object)
      annotation <- 'annotation'
    }
    order_re <- order_receptorgene(object,annotation,allreceptorgenes,allcells,...)
    dat$receptorgene <- factor(dat$receptorgene, levels = order_re)
    dat <- dat[order(dat$receptorgene),]
    dat$interacting_pair <- as.vector(dat$interacting_pair)
    order_pairs <- unique(dat$interacting_pair)
    dat$interacting_pair <- factor(dat$interacting_pair, levels = order_pairs)
  }
  if (!is.null(order.genepair)) {
    dat$interacting_pair <- factor(dat$interacting_pair, levels = order.genepair)
  }
  # plot theme
  mytheme <- theme_bw() + theme(plot.title = element_blank(),#element_text(family="ArialMT", size=25, hjust = 0.5),
                                panel.grid= element_blank(),
                                axis.title = element_blank(),
                                axis.ticks.length = unit(0.3, "cm"),
                                axis.text=element_text(family="ArialMT", size = size.axis.text),
                                axis.text.x=element_text(angle=90, vjust=0.5, hjust=1,size=size.axis.text.x),
                                legend.title=element_text(family="ArialMT", size=size.legend.title),
                                legend.text=element_text(family="ArialMT", size=size.legend.text),
                                legend.background=element_rect(fill="transparent"),
                                plot.margin = unit(c(1,1,1,5), "cm"))
  mycol <- brewer.pal(11,"RdYlBu")
  ##seperate data to avoid  cons memory exhausted.
  if (length(unique(dat$variable)) > max.cell) {
    dat1 <- dat[dat$variable %in% unique(dat$variable)[1:max.cell],]
    dat2 <- dat[dat$variable %in% unique(dat$variable)[max.cell:length(unique(dat$variable))],]
    p1 <- ggplot(dat1, aes(x=variable, y=interacting_pair, colour = means, size=pval_group)) + geom_point() +
      scale_color_gradientn(name='Interaction\nScore', colors=rev(mycol)) +
      scale_size_manual(values = scale.pt.size) +
      coord_flip()+
      mytheme
    p2 <- ggplot(dat2, aes(x=variable, y=interacting_pair, colour = means, size=pval_group)) + geom_point() +
      scale_color_gradientn(name='Interaction\nScore', colors=rev(mycol)) +
      scale_size_manual(values = scale.pt.size) +
      coord_flip()+
      mytheme
    p <- p1+p2
    #plist[[li]] <- p
    if (!is.null(save.pdf)) {
      width1 <- 60 * length(unique(dat1$interacting_pair))/350 +20
      height1 <- ceiling(max.cell/5) + 6
      pdf(save.pdf,width1,height1)
      print(p1)
      print(p2)
      dev.off()
    }
  }else{
    p <- ggplot(dat, aes(x=variable, y=interacting_pair, colour = means, size=pval_group)) + geom_point() +
      scale_color_gradientn(name='Interaction\nScore', colors=rev(mycol)) +
      scale_size_manual(values = scale.pt.size) +
      coord_flip()+
      mytheme
    #plist[[li]] <- p
    if (!is.null(save.pdf)) {
      width <- 60 * length(unique(dat$interacting_pair))/350 +20
      repnum <- length(unique(re))
      height <- ceiling(repnum/5) + 6
      pdf(save.pdf, width, height)
      p <- ggplot(dat, aes(x=variable, y=interacting_pair, colour = means, size=pval_group)) + geom_point() +
        scale_color_gradientn(name='Interaction\nScore', colors=rev(mycol)) +
        coord_flip()+
        mytheme
      print(p)
      dev.off()
    }
  }
  #}
  #return data
  if (return.data) {
    return(dat)
  }else{
    return(p)
  }
}

#' DotPlot modified from Seurat::DotPlot
#' 
#' @export
DotPlot2 <- function (object, features, assay = NULL, col.min = -2.5, col.max = 2.5, 
                      dot.min = 0, dot.scale = 6, color = scPalette_heatmap_bar(), 
                      group.by = NULL, split.by = NULL, scale.by = "radius", scale.min = NA, 
                      scale.max = NA, y.size = 12, x.size = 12, x.vjust = 0.5, 
                      scale = T, x.hjust = 1, x.angle = 90, legend.color.name = "Average\nExpression", 
                      cols = c("lightgrey", "blue")) 
{
  if (!is.null(assay)) {
    DefaultAssay(object = object) <- assay
  }
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  data.features <- FetchData(object = object, vars = features)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)
  }
  else {
    object[[group.by, drop = TRUE]]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]]
    if (length(x = unique(x = splits)) > length(x = cols)) {
      stop("Not enought colors for the number of groups")
    }
    cols <- cols[1:length(x = unique(x = splits))]
    names(x = cols) <- unique(x = splits)
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
                        "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, 
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == 
                                                     x, "avg.exp"]
                             if (scale) {
                               data.use <- scale(x = data.use)
                               data.use <- MinMax(data = data.use, min = col.min, 
                                                  max = col.max)
                             }
                             else {
                               data.use <- log(x = data.use)
                             }
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (!is.null(x = split.by)) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
                                         breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = rev(x = features))
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (!is.null(x = split.by)) {
    splits.use <- vapply(X = strsplit(x = as.character(x = data.plot$id), 
                                      split = "_"), FUN = "[[", FUN.VALUE = character(length = 1L), 
                         2)
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = is.null(x = split.by), yes = "avg.exp.scaled", 
                     no = "colors")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", 
                                                        y = "id")) + geom_point(mapping = aes_string(size = "pct.exp", 
                                                                                                     color = color.by)) + geom_point(mapping = aes_string(size = "pct.exp"), 
                                                                                                                                     color = "black", pch = 21) + scale.func(range = c(0, 
                                                                                                                                                                                       dot.scale), limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(), 
                                                                                                                                                                                                                                             axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) + 
    labs(x = "Features", y = ifelse(test = is.null(x = split.by), 
                                    yes = "Identity", no = "Split Identity")) + cowplot::theme_cowplot()
  plot <- plot + scale_color_gradientn(name = legend.color.name, 
                                       colors = rev(color)) + theme_bw() + xlab("") + ylab("") + 
    theme(axis.text.x = element_text(size = x.size, color = "black", 
                                     vjust = x.vjust, hjust = x.hjust, angle = x.angle), 
          axis.text.y = element_text(size = y.size, color = "black"), 
          panel.grid.major = element_line(size = 0.1, linetype = "dotted", 
                                          colour = "grey"))
  return(plot)
}

