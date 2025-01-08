files <- list.files('./data-raw/')

txtfiles <- files[grepl('txt$',files)]
orthologslist <- lapply(txtfiles, function(x) {
  df <- read.table(paste0('./data-raw/', x), sep = '\t')
})
versions <- gsub('^.*version_|.txt$','', txtfiles)
names(orthologslist) <- paste0('v',versions)
##save rda
# usethis::use_data(orthologslist,
#                   overwrite = TRUE,
#                   internal = TRUE)
