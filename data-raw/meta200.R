## meta data for example
meta200 <- read.table(file = "./data-raw/meta200.txt")
##save rda
usethis::use_data(meta200, overwrite = T)
