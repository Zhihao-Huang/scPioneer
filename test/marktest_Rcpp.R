
###test different cell number
dim(pbmc)
#[1] 13714  2638
data1 <- bench::mark(test1 <- angrycell(pbmc, test.method = 1),
            test2 <- angrycell(pbmc, test.method = 2),
            test3 <- angrycell(pbmc, test.method = 3),
            iterations = 30,check = F)
nasal <- readRDS('/BGFS1/home/huangzh/workspace/project/nasal_polyp/result/10.result_without_NP7/090.merge_cluster_result/SeuratS4.rds')
Idents(nasal) <- nasal$Annotation
head(Idents(nasal))
dim(nasal)
#[1]  28764 122664
data2 <- bench::mark(test1 <- angrycell(nasal, test.method = 1),
            test2 <- angrycell(nasal, test.method = 2),
            test3 <- angrycell(nasal, test.method = 3),
            iterations = 30,check = F)
ov <- readRDS('/BGFS1/home/huangzh/workspace/project/ovarian_cancer/result/RNA/160.summarize_clustering/SeuratS4.rds')
Idents(ov) <- ov$Annotation
DimPlot(ov, group.by = 'Annotation')
dim(ov)
#[1]  27157 237310
data3 <- bench::mark(test1 <- angrycell(ov, test.method = 1),
            test2 <- angrycell(ov, test.method = 2),
            test3 <- angrycell(ov, test.method = 3),
            iterations = 30,check = F)


###test different lib
data4 <- bench::mark(test1 <- angrycell_v1.1.7(pbmc, test.method = 1,select.db = c('db_self','db_lung_zzm','db_immune_qss','db_cellmarker',
                                                                            'db_ovarian_cancer_T','lung_markers_zhangzemin_Cell2019')),
                     test2 <- angrycell_v1.1.7(pbmc, test.method = 2,select.db = c('db_self','db_lung_zzm','db_immune_qss','db_cellmarker',
                                                                            'db_ovarian_cancer_T','lung_markers_zhangzemin_Cell2019')),
                     test3 <- angrycell_v1.1.7(pbmc, test.method = 3,select.db = c('db_self','db_lung_zzm','db_immune_qss','db_cellmarker',
                                                                            'db_ovarian_cancer_T','lung_markers_zhangzemin_Cell2019')),
                     iterations = 3,check = F)
data5 <- bench::mark(test1 <- angrycell_v1.1.7(pbmc, test.method = 1,select.db = c('db_self')),
                     test2 <- angrycell_v1.1.7(pbmc, test.method = 2,select.db = c('db_self')),
                     test3 <- angrycell_v1.1.7(pbmc, test.method = 3,select.db = c('db_self')),
                     iterations = 3,check = F)
