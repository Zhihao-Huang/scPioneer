#' Marker data frame to input data for add.db
#' 
#' @examples 
#' mdf <- data.frame(cellName = c('Naive_CD4T','Naive_CD4T','Naive_CD4T','Epithelial','Epithelial'),
#' geneSymbol = c('CD3D','CD4','CCR7','EPCAM','KRT8'), stringsAsFactors = F)
#' user_db <- add_db_df(mdf)
#' 
#' @export
add_db_df <- function (mdf = NULL) {
  df <- mdf %>% group_by(cellName) %>% 
    summarise(geneSymbol = paste0(geneSymbol, collapse = ","))
  return(df)
}

#' Marker list to input data for add.db
#' 
#' @examples 
#' mlist <- list('Naive CD4 T cell' = c('CD3D','CD4','CCR7'),'Epithelial' = c('EPCAM','KRT8'))
#' user_db <- add_db_list(mlist)
#' 
#' @export
add_db_list <- function (mlist = NULL) {
  mlist_str <- list()
  for (n in names(mlist)) {
    str <- paste0(mlist[[n]], collapse = ",")
    mlist_str[[n]] <- str
  }
  addm <- as.data.frame(do.call(rbind, mlist_str), stringsAsFactors = F)
  addm$Celltype <- rownames(addm)
  colnames(addm) <- c("geneSymbol", "cellName")
  addm <- addm[, c("cellName", "geneSymbol")]
  return(addm)
}

#' Search markers or cell types.
#' @export
angrycelldb <- function (genename = NULL, cellname = NULL, use.specific.marker = F) {
  mlist <- c(mlist1,mlist2,mlist3)
  if (!is.null(cellname)) {
    pos <- cellname %in% names(mlist)
    if (sum(pos) > 0) {
      cells <- mlist[cellname[pos]]
      cells <- lapply(cells, function(x) paste(x, collapse = ', '))
      refs <- Rlist[cellname[pos]]
      cellmat <- do.call(rbind, cells)
      refmat <- do.call(rbind, refs)
      mat <- cbind(cellmat,refmat)
      colnames(mat) <- c('Marker gene','Referece')
    }else{
      similar_type <- agrep(cellname, names(mlist), max = 1, ignore.case = T, value = T)
      stop(paste0('All the inputed cell types are not in AngryCellDB. Do you want to search: \n',
                  paste(similar_type, collapse = ', '), '?'))
    }
    return(mat)
  }else if (!is.null(genename)) {
    if (use.specific.marker) {
      genelist <- vector("list", length = length(genename))
      names(genelist) <- genename
      for (g in genename){
        for (i in names(slist)) {
          if (g %in% slist[[i]]) {
            genelist[[g]] <- c(genelist[[g]], i)
          }
        }
      }
      genelist <- lapply(genelist, function(x) paste(x, collapse = ', '))
      genemat <- do.call(rbind, genelist)
      colnames(genemat) <- c('Cell type')
    }else{
      genelist <- vector("list", length = length(genename))
      names(genelist) <- genename
      for (g in genename){
        for (i in names(mlist)) {
          if (g %in% slist[[i]]) {
            genelist[[g]] <- c(genelist[[g]], i)
          }
        }
      }
      genelist <- lapply(genelist, function(x) paste(x, collapse = ', '))
      genemat <- do.call(rbind, genelist)
      colnames(genemat) <- c('Cell type')
    }
    return(genemat)
  }else{
    message('Please input cellname or genename for searching.')
  }
}

#' Show all cell types in tree。
#'
#' @export

List_all_celltypes <- function(outfile = NULL, type = "text"){
  input <- c()
  for (i in names(mindrlist1)) {
    type1 <- paste0('# ',i)
    input <- c(input, type1)
    if (length(mindrlist1[[i]]) > 0) {
      type2list <-mindrlist1[[i]]
      for (j in type2list) {
        input <- c(input,  paste0('## ', j))
        if (length(mindrlist2[[j]]) > 0) {
          input <- c(input,  paste0('### ', mindrlist2[[j]]))
        }
      }
    }
  }
  mindr::mm(from = input, to = outfile, type = type, root = "Cell types")
}

#' seurat object
#'
#' a seurat object for example
#'
#' @examples
#' pbmc <- data('pbmc')
"pbmc"

#' meta data
#'
#' @examples
#' meta200 <- data('meta200')
"meta200"

#' Template of uploading markers for AngryCellDB v1.2.0
#' @examples
#' Template()
#' @export
Template <- function(version = 'v1.2.0') {
  if (version == 'v1.2.0') {
    h <- head(templv1.2.0,3)
  }
  return(h)
}
#' Upload markers for AngryCellDB v1.2.0
#' @examples
#' upload_marker('/BGFS1/projectdata/jasper/database/marker_collection/marker_collection_template_v1.2.0.txt', author = "ShaoYu Song")
#' @export
upload_marker <- function (txtfilename, 
          author = "anonymous", 
          outdir = "/BGFS1/projectdata/jasper/database/marker_collection", 
    version = NULL) 
{
    txtfile <- read.csv(sep = "\t", check.names = F, stringsAsFactors = F, 
        file = txtfilename)
    d <- date()
    user_links <- user_links
    user_Rlist <- user_Rlist
    user_slist <- user_slist
    user_mindrlist1 <- user_mindrlist1
    user_mindrlist2 <- user_mindrlist2
    user_mlist1 <- user_mlist1
    user_mlist2 <- user_mlist2
    user_mlist3 <- user_mlist3
    user_mlist <- c(user_mlist1, user_mlist2, user_mlist3)
    if (is.null(version)) {
        version <- packageVersion("angrycell")[[1]]
    }
    sharefile <- paste0(outdir, "/marker_collection_v", version, 
        ".R")
    message(paste0("Adding marker genes to file: \n", sharefile))
    if (!file.exists(sharefile)) {
        file.create(sharefile)
    }
    source(sharefile)
    write(paste0("\n\n\n# Upload date: ", d, "; Author: ", author, 
        "."), file = sharefile, append = TRUE)
    if (sum(txtfile$species == "Homosapiens") != 0) {
        Human_marker <- txtfile[txtfile$species == "Homosapiens", 
            ]
        for (i in 1:nrow(Human_marker)) {
            refpos <- sapply(user_links, function(x) Human_marker[i, 
                "references"] %in% x)
            if (sum(refpos) > 0) {
                refname <- paste0("lks", which(refpos)[1])
            }
            else {
                numlks <- length(names(user_links)) + 1
                refname <- paste0("lks", numlks)
                user_links[[refname]] <- Human_marker[i, "references"]
            }
            sametype <- Human_marker[i, "celltypes"] %in% c(names(user_mlist))
            if (!sametype) {
                user_Rlist[[Human_marker[i, "celltypes"]]] <- Human_marker[i, 
                  "references"]
            }
            sameref <- Human_marker[i, "references"] == user_Rlist[[Human_marker[i, 
                "celltypes"]]]
            if (sametype & !sameref) {
                Human_marker[i, "celltypes"] <- paste0(Human_marker[i, 
                  "celltypes"], "_", refname)
                sametype <- Human_marker[i, "celltypes"] %in% 
                  c(names(user_mlist))
                if (sametype) {
                  message(paste0("Existed Cell type (", i, "): ", 
                    Human_marker[i, "celltypes"], ", pass."))
                  next
                }
            }
            if (sametype & sameref) {
                message(paste0("Existed Cell type (", i, "): ", 
                  Human_marker[i, "celltypes"], ", pass."))
                next
            }
            user_mlist[[Human_marker[i, "celltypes"]]] <- 1
            user_Rlist[[Human_marker[i, "celltypes"]]] <- user_links[[refname]]
            message(paste0("Recording cell type (", i, "): ", 
                Human_marker[i, "celltypes"]))
            write(paste0("user_mlist", Human_marker[i, "layers"], 
                "[['", Human_marker[i, "celltypes"], "']] <- c('", 
                gsub(",", "','", Human_marker[i, "markers"]), 
                "')"), file = sharefile, append = TRUE)
            write(paste0("user_slist[['", Human_marker[i, "celltypes"], 
                "']] <- c('", gsub(",", "','", Human_marker[i, 
                  "specific_markers"]), "')"), file = sharefile, 
                append = TRUE)
            write(paste0("user_links[['", refname, "']] <- '", 
                Human_marker[i, "references"], "'"), file = sharefile, 
                append = TRUE)
            write(paste0("user_Rlist[['", Human_marker[i, "celltypes"], 
                "']] <- user_links[['", refname, "']]"), file = sharefile, 
                append = TRUE)
            pos2 <- Human_marker[i, "up_layers"] %in% names(user_mindrlist1)
            pos3 <- Human_marker[i, "up_layers"] %in% names(user_mindrlist2)
            if (pos2) {
                user_mindrnum <- 1
                write(paste0("user_mindrlist", user_mindrnum, 
                  "[['", Human_marker[i, "up_layers"], "']] <- c(user_mindrlist", 
                  user_mindrnum, "[['", Human_marker[i, "up_layers"], 
                  "']],'", Human_marker[i, "celltypes"], "')"), 
                  file = sharefile, append = TRUE)
            }
            else if (pos3) {
                user_mindrnum <- 2
                write(paste0("user_mindrlist", user_mindrnum, 
                  "[['", Human_marker[i, "up_layers"], "']] <- c(user_mindrlist", 
                  user_mindrnum, "[['", Human_marker[i, "up_layers"], 
                  "']],'", Human_marker[i, "celltypes"], "')"), 
                  file = sharefile, append = TRUE)
            }
            else if (!pos2 & (!pos3)) {
                user_mindrnum <- 1
                write(paste0("user_mindrlist", user_mindrnum, 
                  "[['", Human_marker[i, "celltypes"], "']] <- c()"), 
                  file = sharefile, append = TRUE)
            }
        }
    }
    if (sum(txtfile$species == "mm10") != 0) {
        Mouse_marker <- txtfile[txtfile$species == "mm10", ]
        for (i in 1:nrow(Mouse_marker)) {
            write(paste0("mm10list", "[['", Mouse_marker[i, "celltypes"], 
                "']] <- c('", gsub(",", "','", Mouse_marker[i, 
                  "markers"]), "')"), file = sharefile, append = TRUE)
        }
    }
    write("user_mlist <- c(user_mlist1, user_mlist2, user_mlist3)", 
        file = sharefile, append = TRUE)
    write("# Upload accomplished.\n\n\n", file = sharefile, append = TRUE)
    message("Data uploaded complete.")
}

#' upload marker
#' 
#' @param cellname cell type name
#' @param markers general expressing genes.
#' @param Definite_markers specific markers of this cell type.
#' @param author user who uploads this data.
#' 
#' @examples
#' uploaddb('CD4+ naive T cell', markers = c('CD4','LEF1','TCF7','CCR7','SELL'),
#' Definite_markers = c('CD4','CCR7'), author = 'Zhihao Huang')
#' 
#' @export
uploaddb <- function (cellname, markers = "_", Definite_markers = "_", author = "User") 
{
  version <- packageVersion("angrycell")[[1]]
  sharefile <- paste0("/BGFS1/projectdata/jasper/database/marker_collection/marker_collection_", 
                      version, ".txt")
  if (!file.exists(sharefile)) {
    file.create(sharefile)
  }
  if (length(markers) > 1) {
    mk <- paste0(markers, collapse = ",")
  }
  if (length(Definite_markers) > 1) {
    mkd <- paste0(Definite_markers, collapse = ",")
  }
  time <- format(Sys.time(), "%a %b %d %X %Y")
  addm <- t(as.matrix(c(cellname, mk, mkd, author, time)))
  db <- read.csv(file = sharefile, sep = "\t", header = T)
  colnames(addm) <- colnames(db)
  db <- as.data.frame(rbind(db, addm), stringsAsFactors = F)
  write.table(db, file = sharefile, quote = F, sep = "\t", 
              row.names = F)
  print("Data uploaded complete.")
}