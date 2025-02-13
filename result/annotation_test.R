### Create a samplelist including samplename and datadir.
df <- data.frame(samplename = 'PBMC', datadir = './data-raw/filtered_gene_bc_matrices/hg19/')
write.table(df, file = './result/samplelist.txt', sep = '\t')


## load data
library(scPioneer)
source('./R/annotation.R')
### set arguments
param <- PHASE1_run_Seurat_v5_QC_clustering_param_template()
param$samplelist <- './result/samplelist.txt'
param$outdir <- './result/'
param$min.features <- 200
param$max.mt <- 10
param$min.cells <- 3
param$species <- 'Human'
param$normalize_method <- 'SCT'
param$detect.doublet <- "scDblFinder"
param$return.plot <- T

### run
resultlist <- PHASE1_run_Seurat_v5_QC_clustering(param)

names(resultlist)
# object plotlist
patchwork::wrap_plots(resultlist$plotlist)

saveRDS(resultlist$object, file = './result/pbmc_human.rds')
# get DEGs by FindAllMarkers
Idents(pbmc) <- pbmc$seurat_clusters
markers <- FindAllMarkers(pbmc, logfc.threshold = 0.5, test.use = 'MAST', only.pos = T)
write.table(markers, sep = '\t', file = './result/deg_human.txt')

### Load preprecessed object from PHASE 1.
pbmc <- resultlist$object

### step 2
library(Seurat)
library(dplyr)
library(ggplot2)
source('./R/annotation.R')
source('./R/DimPlot_idx.R')
pbmc <- readRDS('./result/pbmc_human.rds')
plist <- list()
### Perform annotation by reference (SingleR)
obj <- annocell(pbmc_small, species = 'Human', method = 'SingleR', assay = 'SCT')
p1 <- DimPlot_idx(obj) + ggtitle('singleR')
ggsave('./result/anno_singler.png', device = 'png', plot = p1, width = 5.2,height = 3.25)
plist[['singler']] <- p1
### Perform annotation by LLM model 
markers <- read.table('./result/deg_human.txt')
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
# OpenAI GPT-3.5-turbo
obj <- annocell(pbmc, species = 'Human', 
                method = 'llm', 
                llm_function = 'openai',
                model = "gpt-3.5-turbo",
                DE = top10,  seed = 1234, sep = ': ', as.order = F,
                base_url = "http://chatapi.littlewheat.com/v1",
                api_key = 'sk-7OpB4wXOH9fmCngPJAg0MPeqIU0JYxvVnSvk10pzZX0Bj1H6')
p2 <- DimPlot_idx(obj) + ggtitle('gpt-3.5-turbo')
ggsave('./result/anno_gpt-3.5-turbo.png', device = 'png', plot = p2, width = 6.25,height = 3.25)
plist[['gpt']] <- p2

# DeepSeek R1
obj <- annocell(pbmc, species = 'Human', method = 'llm', model = "deepseek-reasoner",
                DE = top10,  seed = 1234, sep = ': ', as.order = F,
                base_url = 'https://api.deepseek.com/v1',
                api_key = 'sk-dfa7e4495c174050a82d81e6f82a9d14')
p3 <- DimPlot_idx(obj) + ggtitle('deepseek-chat')
ggsave('./result/anno_DeepSeekR1.png', device = 'png', plot = p3, width = 7,height = 3.25)
plist[['DeepSeekR1']] <- p3

text <- anno_ellmer(top20, method = 'llm', model = "deepseek-reasoner",
                    seed = 1234, sep = '\\. ', as.order = F,
                    base_url = 'https://api.deepseek.com/v1',
                    api_key = 'sk-dfa7e4495c174050a82d81e6f82a9d14')
write.table(annodf20, sep = '\t', quote = F, 
            file = paste0(outdir, 'annodf20.txt'))


### Perform annotation by top markers
markerdf <- data.frame(celltypes = c('T','T','NK','NK','Mono', 'Mono','DC','DC','B','B','Platelet'), 
                       markers = c('CD3D','CD3E','NCAM1','NKG7','CD14','FCGR3A','CST3','CD1C','CD79A','MS4A1','PPBP'))
colnames(markerdf)
# "celltypes" "markers"
obj <- annocell(pbmc, species = 'Human', method = 'topgene', markerdf = markerdf)
p4 <- DimPlot_idx(obj) + ggtitle('top gene')
ggsave('./result/anno_topgene.png', device = 'png', plot = p4, width = 5,height = 3.25)
plist[['topgene']] <- p4
saveRDS(plist, file= './result/plist.rds')

# angrycell + Deepseek
obj <- annocell(pbmc, species = 'Human', method = 'angrycell', db = 'openai', assay = 'SCT',
                DE = top10,model = "gpt-3.5-turbo", seed = 1234, 
                sep = ': ', as.order = F,
                base_url = "http://chatapi.littlewheat.com/v1",
                api_key = 'sk-7OpB4wXOH9fmCngPJAg0MPeqIU0JYxvVnSvk10pzZX0Bj1H6')
p3 <- DimPlot_idx(obj) + ggtitle('angrycell + deepseek-chat')
p3

sessionInfo()


