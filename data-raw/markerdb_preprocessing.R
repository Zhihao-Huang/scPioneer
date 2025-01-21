#### panglaodb
library(dplyr)
## NOT RUN
if (F) {
## prepare additional database
panglaodb <- read.table(sep = '\t', header = T,
                        file = './data-raw/PanglaoDB_markers_27_Mar_2020.tsv.gz')

# mm10
orthologs <- biospiper::Get_orthologs_mouse_human(version = 98,
                                                  remove.duplicated = F,
                                                  only.one.to.one = F,
                                                  using.local.file = T)
orthologs <- orthologs[!duplicated(orthologs$HGNC_symbol) & orthologs$HGNC_symbo != '',]
panglaodb_mm <- panglaodb[grepl('Mm',panglaodb$species),]
colnames(panglaodb_mm)[2] <- 'HGNC_symbol'
panglaodb_mm <- left_join(panglaodb_mm, orthologs)

dim(panglaodb_mm)
#[1] 6087   17
sum(is.na(panglaodb_mm$MGI_symbol))
#[1] 524
panglaodb_mm <- panglaodb_mm[!is.na(panglaodb_mm$MGI_symbol),]
# only specific markers; 
panglaodb_mm <- panglaodb_mm[!duplicated(panglaodb_mm$MGI_symbol),]
panglaodb_mm <- panglaodb_mm %>% group_by(cell.type) %>% 
  summarise(geneSymbol = paste(MGI_symbol, collapse = ','))
colnames(panglaodb_mm)[1] <- 'cellName'
panglaodb_mm <- as.data.frame(panglaodb_mm)
rownames(panglaodb_mm) <- panglaodb_mm$cellName
#write.table(panglaodb_mm, quote = F, sep = '\t', 
#            file = './data-raw/PanglaoDB_markers_27_Mar_2020_mm.txt')

# hs
panglaodb_hs <- panglaodb[grepl('Hs',panglaodb$species),]
colnames(panglaodb_hs)[2] <- 'HGNC_symbol'
dim(panglaodb_hs)
#[1] 6034   17
sum(is.na(panglaodb_hs$HGNC_symbol))
#[1] 0
# only specific markers; 
panglaodb_hs <- panglaodb_hs[!duplicated(panglaodb_hs$HGNC_symbol),]
panglaodb_hs <- panglaodb_hs %>% group_by(cell.type) %>% 
  summarise(geneSymbol = paste(HGNC_symbol, collapse = ','))
colnames(panglaodb_hs)[1] <- 'cellName'
panglaodb_hs <- as.data.frame(panglaodb_hs)
rownames(panglaodb_hs) <- panglaodb_hs$cellName
#write.table(panglaodb_hs, quote = F, sep = '\t', 
#            file = './data-raw/PanglaoDB_markers_27_Mar_2020_hs.txt')


  # databases from angrycell
  lazyLoad('/BGFS1/usr/share/packages/R/Rlibs-3.6.1/angrycell/R/sysdata')
  sysdata <- ls()
  paste(sysdata,collapse = ', ')
  angrycelldata <- list(A_M, addm_qss, D_M, db_self, db_self_mm10, db_user_mm10,
                        db_users, db1, db2, db3, db4, dblayers, dbref, dbspf, HKgene,
                        homolog_Human_mouse, links, mindrlist1, mindrlist2, mlist1,
                        mlist2, mlist3, MSigDBlist, MSigDBnames, Rlist, slist, 
                        templv1.2.0, user_links, user_mindrlist1, user_mindrlist2, 
                        user_mlist1, user_mlist2, user_mlist3, user_Rlist, user_slist,
                        panglaodb_hs,panglaodb_mm)
  names(angrycelldata) <- c('A_M','addm_qss','D_M','db_self','db_self_mm10','db_user_mm10',
                            'db_users','db1','db2','db3','db4','dblayers','dbref','dbspf','HKgene',
                            'homolog_Human_mouse','links','mindrlist1','mindrlist2','mlist1',
                            'mlist2','mlist3','MSigDBlist','MSigDBnames','Rlist','slist',
                            'templv1.2.0','user_links','user_mindrlist1','user_mindrlist2',
                            'user_mlist1','user_mlist2','user_mlist3','user_Rlist','user_slist',
                            'panglaodb_hs','panglaodb_mm')
  saveRDS(angrycelldata, file = '/BGFS1/projectdata/jasper/project/pipeline_abio_huangzhihao/scPioneer/scpioneer_v1.0.0/scPioneer/data-raw/angrycelldb.rds')
  
}
angrycelldata <- readRDS('./data-raw/angrycelldb.rds')
print(names(angrycelldata))
cat(paste(paste0(names(angrycelldata), ' <- angrycelldata$', names(angrycelldata)),collapse = '\n'))
A_M <- angrycelldata$A_M
addm_qss <- angrycelldata$addm_qss
D_M <- angrycelldata$D_M
db_self <- angrycelldata$db_self
db_self_mm10 <- angrycelldata$db_self_mm10
db_user_mm10 <- angrycelldata$db_user_mm10
db_users <- angrycelldata$db_users
db1 <- angrycelldata$db1
db2 <- angrycelldata$db2
db3 <- angrycelldata$db3
db4 <- angrycelldata$db4
dblayers <- angrycelldata$dblayers
dbref <- angrycelldata$dbref
dbspf <- angrycelldata$dbspf
HKgene <- angrycelldata$HKgene
homolog_Human_mouse <- angrycelldata$homolog_Human_mouse
links <- angrycelldata$links
mindrlist1 <- angrycelldata$mindrlist1
mindrlist2 <- angrycelldata$mindrlist2
mlist1 <- angrycelldata$mlist1
mlist2 <- angrycelldata$mlist2
mlist3 <- angrycelldata$mlist3
MSigDBlist <- angrycelldata$MSigDBlist
MSigDBnames <- angrycelldata$MSigDBnames
Rlist <- angrycelldata$Rlist
slist <- angrycelldata$slist
templv1.2.0 <- angrycelldata$templv1.2.0
user_links <- angrycelldata$user_links
user_mindrlist1 <- angrycelldata$user_mindrlist1
user_mindrlist2 <- angrycelldata$user_mindrlist2
user_mlist1 <- angrycelldata$user_mlist1
user_mlist2 <- angrycelldata$user_mlist2
user_mlist3 <- angrycelldata$user_mlist3
user_Rlist <- angrycelldata$user_Rlist
user_slist <- angrycelldata$user_slist
panglaodb_hs <- angrycelldata$panglaodb_hs
panglaodb_mm <- angrycelldata$panglaodb_mm
# orthologs genesets
files <- list.files('./data-raw/')
txtfiles <- files[grepl('^orthologs_mouse',files)]
orthologslist <- lapply(txtfiles, function(x) {
  df <- read.table(paste0('./data-raw/', x), sep = '\t')
})
versions <- gsub('^.*version_|.txt$','', txtfiles)
names(orthologslist) <- paste0('v',versions)

# symbol ensembl genelist
files <- list.files('./data-raw/sym_ebl/Human/')
txtfiles <- files[grepl('txt$',files)]
humanlist <- lapply(txtfiles, function(x) {
  df <- read.table(paste0('./data-raw/sym_ebl/Human/', x), sep = '\t')
})
versions <- gsub('^.*version_|.txt$','', txtfiles)
names(humanlist) <- paste0('v',versions)
files <- list.files('./data-raw/sym_ebl/Mouse/')
txtfiles <- files[grepl('txt$',files)]
mouselist <- lapply(txtfiles, function(x) {
  df <- read.table(paste0('./data-raw/sym_ebl/Mouse/', x), sep = '\t')
})
versions <- gsub('^.*version_|.txt$','', txtfiles)
names(mouselist) <- paste0('v',versions)
symbol_ensembl_list <- list(Human = humanlist, Mouse = mouselist)

##save rda
usethis::use_data(A_M, addm_qss, D_M, db_self, db_self_mm10, db_user_mm10,
                  db_users, db1, db2, db3, db4, dblayers, dbref, dbspf, HKgene,
                  homolog_Human_mouse, links, mindrlist1, mindrlist2, mlist1,
                  mlist2, mlist3, MSigDBlist, MSigDBnames, Rlist, slist, 
                  templv1.2.0, user_links, user_mindrlist1, user_mindrlist2, 
                  user_mlist1, user_mlist2, user_mlist3, user_Rlist, user_slist,
                  panglaodb_hs,panglaodb_mm, 
                  overwrite = TRUE,
                  internal = TRUE)


