### Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
source('./R/annotation.R')
source('./R/DimPlot_idx.R')
# Load data
pbmc <- readRDS('./result/pbmc_human.rds')
# perform deg analysis
Idents(pbmc) <- pbmc$seurat_clusters
deg <- FindAllMarkers(pbmc, test.use = 'MAST', only.pos = T)
write.table(deg, sep = '\t', quote = F, 
            file = './result/deg_human.txt')

### Perform annotation by LLM model 
markers <- read.table('./result/deg_human.txt')

top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)

signgene <- markers %>% group_by(cluster) %>% filter(p_val_adj < 0.05) 
# T cell
unique(signgene$cluster[signgene$gene == 'CD3D'])
# [1] 0 2 4 7 8
unique(signgene$cluster[signgene$gene == 'CD8B'])
# [1] 4 7 8
unique(signgene$cluster[signgene$gene == 'GZMK'])
# [1] 4 7
unique(signgene$cluster[signgene$gene == 'IL7R'])
# [1] 0 2 4
unique(signgene$cluster[signgene$gene == 'IL2RA'])
# [1] 2
unique(signgene$cluster[signgene$gene %in% c('CCR7')])
# [1] 0 8
# NK
unique(signgene$cluster[signgene$gene %in% c('NCR3')])
# [1] 4 6
unique(signgene$cluster[signgene$gene %in% c('NKG7')])
# [1] 4 6 7
# Monocyte
unique(signgene$cluster[signgene$gene %in% c('CD14')])
# [1] 1
# Platelets
unique(signgene$cluster[signgene$gene %in% c('PPBP')])
# [1] 5
# B cell
unique(signgene$cluster[signgene$gene %in% c('CD79B')])
# [1] 3
unique(signgene$cluster[signgene$gene %in% c('CD1C')])
# [1] 3 9
unique(signgene$cluster[signgene$gene %in% c('FCER1A')])
# [1] 9

Manual_anno <- c('Naive CD4+ T cells', 'Monocytes', 'Memory CD4+ T cells',
                 'B cells', 'Effector memory CD8+ T cells','Platelets',
                 'NK cells', 'Cytotoxic CD8+ T cells','CD8+ T cells',
                 'Dendritic cells')

annolist <- list()
#### API-based
# OpenAI GPT-3.5-turbo
annolist[['GPT-3.5-turbo']] <-  anno_ellmer(top10, tissuename = 'PBMC',
                                            llm_function = 'openai', 
                                            model = "gpt-3.5-turbo",
                                            seed = 1234, sep = ': ', as.order = F,
                                            base_url = 'http://chatapi.littlewheat.com/v1',
                                            api_key = 'sk-7OpB4wXOH9fmCngPJAg0MPeqIU0JYxvVnSvk10pzZX0Bj1H6')

# DeepSeek V3 (671b)
annolist[['DeepSeekV3']] <- anno_ellmer(top10, tissuename = 'PBMC',
                                        llm_function = 'deepseek', 
                                        model = "deepseek-chat",
                                        seed = 1234, sep = ': ', as.order = F,
                                        base_url = 'https://api.deepseek.com/v1',
                                        api_key = 'sk-dfa7e4495c174050a82d81e6f82a9d14')


# DeepSeek R1 (671b)
annolist[['deepseek-r1:671b']] <- anno_ellmer(top10, tissuename = 'PBMC',
                                        llm_function = 'deepseek', 
                                        model = "deepseek-reasoner",
                                        #return.answer = T,
                                        seed = 1234, sep = ': ', as.order = F,
                                        base_url = 'https://api.deepseek.com/v1',
                                        api_key = 'sk-dfa7e4495c174050a82d81e6f82a9d14')


#### local models (deepseek-r1:70b/7b, others: 7-8b)
#### ./ollama serve # activate ollama models 
models <- c('deepseek-r1:70b','deepseek-r1:7b',
            'llama3.1','qwen2.5',
            'dolphin3',  'gemma2', 
            'mistral', 'falcon3', 'exaone3.5')
# bad models: #'olmo2',#, 'tulu3','llava'
### Res memory and CPUs: 42.4GB + 50 threads, 5-8GB + 50 threads.
# running models
times <- c()
for (i in models) {
  start.time <- Sys.time()
  annolist[[i]] <- anno_ellmer(top10, llm_function = 'ollama', ollama_model = i, 
                               tissuename = 'PBMC', return.answer = F,
                               prompts.add = 'Only return the name of cell type for each row.\n Do not interpret.',
                               rm_str = "^.*think\\>|.*cell types based on marker associations:|\n\n|\n\nEach cluster's composition.*$|The analysis of .* based on marker associations:
\n\n|\\*\\*|Here are the .* marker sets:\n\n|Here are the .* for each row:\n|Here are the identified cell types for each row:\n\n|For each row of cell types .* returned names of cell types:\n",
                               seed = 1234, sep = c('\\: ','\\:', '\\. '), as.order = F)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  times <- c(times, time.taken)
  
}
saveRDS(annolist, file = './result/annodf_list.rds')
names(times) <- models
times
# deepseek-r1:7b  deepseek-r1:70b        llama3.1         qwen2.5        dolphin3
# 3.903418 mins       15.770963 mins      21.426141 secs      23.920832       48.966797
# mistral         falcon3       exaone3.5     gemma2        llava
# 34.219434       37.128616       53.338188  3.69666 mins   18.80901 seconds


#### output results as data frame
annolist <- readRDS('./result/annodf_list.rds')
# ds7b cluster0: monocytes/macrophages
annolist2 <- lapply(names(annolist), function(x) {
  annolist[[x]]$model <- x
  annolist[[x]]
})
names(annolist2) <- names(annolist)
annolist2[['deepseek-r1:7b']] <- rbind(data.frame(Orig_Idents = '0',
                                                  Celltype_predicted = '-',
                                                  model = 'deepseek-r1:7b'),annolist2[['deepseek-r1:7b']])
annolist2[['mistral']]$Orig_Idents <- as.numeric(annolist2[['mistral']]$Orig_Idents) - 1
annolist2[['mistral']] <- rbind(annolist2[['mistral']], 
                                       data.frame(Orig_Idents = 9,
                                                  Celltype_predicted = '-',
                                                  model = 'mistral'))


annodf <- sapply(annolist2, function(x) x$Celltype_predicted)
annodf <- as.data.frame(annodf)
annodf$Manual <- Manual_anno
annodf <- annodf[, c('Manual', 'GPT-3.5-turbo', 'DeepSeekV3',
                     'deepseek-r1:671b',
                     models)]
colnames(annodf)[5:13] <- paste0(colnames(annodf)[5:13],':7b')
write.csv(annodf, file = './result/anno_13_llm_models_local.csv')

sessionInfo()



