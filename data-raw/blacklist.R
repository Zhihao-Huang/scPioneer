## meta data for example
blacklist <- read.csv(file = "./data-raw/excluded_genelist.csv", row.names = 1)
noise_gene <- as.vector(blacklist)
noise_gene <- unlist(blacklist)
noise_gene <- noise_gene[noise_gene != '']
##save rda
usethis::use_data(blacklist, overwrite = T)
usethis::use_data(noise_gene, overwrite = T)
