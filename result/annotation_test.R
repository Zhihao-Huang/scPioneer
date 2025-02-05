### step 2
library(Seurat)
library(dplyr)
library(ggplot2)
source('./R/annotation.R')
source('./R/DimPlot_idx.R')
pbmc <- readRDS('./result/pbmc_mouse.rds')
plist <- list()
### Perform annotation by reference (SingleR)
obj <- annocell(pbmc, species = 'Human', method = 'SingleR', raw_cluster = 'seurat_clusters')
p1 <- DimPlot_idx(obj) + ggtitle('singleR')
ggsave('./result/anno_singler.png', device = 'png', plot = p1, width = 5.2,height = 3.25)
plist[['singler']] <- p1
### Perform annotation by top markers
markerdf <- data.frame(celltypes = c('T','T','NK','NK','Mono', 'Mono','DC','DC','B','B','Platelet'), 
                       markers = c('CD3D','CD3E','NCAM1','NKG7','CD14','FCGR3A','CST3','CD1C','CD79A','MS4A1','PPBP'))
colnames(markerdf)
# "celltypes" "markers"
obj <- annocell(pbmc, species = 'Human', method = 'topgene',
                markerdf = markerdf, raw_cluster = 'seurat_clusters')
p4 <- DimPlot_idx(obj) + ggtitle('top gene')
ggsave('./result/anno_topgene.png', device = 'png', plot = p4, width = 5,height = 3.25)
plist[['topgene']] <- p4
### Perform annotation by LLM model 
markers <- read.table('./result/deg_mouse.txt')
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
# OpenAI GPT-3.5-turbo
obj <- annocell(pbmc, species = 'Human', method = 'llm', 
                llm_function = 'openai', model = "gpt-3.5-turbo",
                raw_cluster = 'seurat_clusters', assay = 'SCT',
                DE = top10,  seed = 1234, sep = ': ', as.order = F,
                base_url = "http://chatapi.littlewheat.com/v1",
                api_key = '')
p2 <- DimPlot_idx(obj) + ggtitle('gpt-3.5-turbo')
ggsave('./result/anno_gpt-3.5-turbo.png', device = 'png', plot = p2, width = 6.25,height = 3.25)
plist[['gpt']] <- p2

# DeepSeek R1
obj <- annocell(pbmc, species = 'Human', method = 'llm',
                llm_function = 'openai',model = "deepseek-chat",
                raw_cluster = 'seurat_clusters',
                DE = top10,  seed = 1234, sep = ': ', as.order = F,
                base_url = 'https://api.deepseek.com/v1',
                api_key = )
p3 <- DimPlot_idx(obj) + ggtitle('deepseek-chat')
ggsave('./result/anno_DeepSeekR1.png', device = 'png', plot = p3, width = 7,height = 3.25)
plist[['DeepSeekR1']] <- p3


# Ollama
obj <- annocell(pbmc, species = 'Human', method = 'llm',
                llm_function = 'ollama', ollama_model = 'llama3.2',
                raw_cluster = 'seurat_clusters',
                DE = top10,  seed = 1234, 
                prompts = 'Identify cell type matching each row. Only return name of cell type. Do not explain.',
                rm_str = c('Here is the list of cell types matching each row:\n\n'),
                sep = '\\. ', as.order = F)
p5 <- DimPlot_idx(obj) + ggtitle('ollama')
ggsave('./result/anno_ollama.png', device = 'png', plot = p5, width = 8,height = 3.35)
plist[['ollama']] <- p5
lay <- rbind(c(1,2,3),
             c(4,5,5))

gridExtra::grid.arrange(p1,p2,p3,p4,p5, layout_matrix = lay)


saveRDS(plist, file= './result/plist_top10.rds')


# archived function
obj <- annocell(pbmc, species = 'Human', method = 'angrycell', db = 'openai',
                DE = top10, raw_cluster = 'seurat_clusters',model = "gpt-3.5-turbo", seed = 1234,
                base_url = "http://chatapi.littlewheat.com/v1", assay = 'SCT',
                api_key = '')
p6 <- DimPlot_idx(obj) + ggtitle('angrycell + deepseek-chat')
p6

sessionInfo()
