##test for library
##data set 1
load('/BGFS1/home/huangzh/workspace/project/annotation_tools_development/pbmc.Rdata')
library(angrycell)
library(Seurat)
pbmc_anno <- angrycell(pbmc, Celltypes.db_self, min.pct = 0.3)
#accuracy 100%!
pbmc_expr <- angrycell(pbmc,basegene = basegene,HK.model.expr=T, area.frac = 0.4)


##data set 2
load("/BGFS1/home/huangzh/workspace/project/nasal_polyp/result/10.result_without_NP7/020.T_cells/01.Tcells_mean.var.plotmeans_0.0125_3disp_0.3_dim10_final.Rdata")
Idents(SeuratS4) <- newI
date()
test <- angrycell(SeuratS4,Celltypes.db_self,min.pct = 0.3)
date()
##30s 17/18
write.table(file='/BGFS1/home/huangzh/workspace/project/annotation_tools_development/test4_T0.3.txt',test,quote = F,sep='\t')
date()
test_T_expr <- angrycell(SeuratS4,Celltypes.db_self,basegene = basegene,HK.model.expr=T, area.frac = 0.4)
date()
write.table(file='/BGFS1/home/huangzh/workspace/project/annotation_tools_development/test4_T_expr.txt',test_T_expr,quote = F,sep='\t')
date()
test_T_drop <- angrycell(SeuratS4,Celltypes.db_self,basegene = basegene,HK.model.drop=T, area.frac = 0.1)
date()
write.table(file='/BGFS1/home/huangzh/workspace/project/annotation_tools_development/test4_T_drop.txt',test_T_drop,quote = F,sep='\t')
date()
test_T_drop2 <- angrycell(SeuratS4,Celltypes.db_self,basegene = basegene,HK.model.drop2=T)
date()
write.table(file='/BGFS1/home/huangzh/workspace/project/annotation_tools_development/test4_T_drop2.txt',test_T_drop2,quote = F,sep='\t')

