#' Calculate diversity of cell types. used Simpson’s index to characterize the composition of cell types.
#' 
#' @examples 
#' ###Calculate diversity for each sample.
#' library(rstatix)
#' library(ggpubr)
#' diver_t_test(meta200, annotation = 'maingroups',sample = 'Sample',
#' group = 'Group', index = 'simpson',return.plot = T)
#' 
#' @export

diver_t_test <- function (meta, annotation = "Annotation", sample = "Sample", 
          group = "Group", index = c("simpson", "invsimpson", "shannon"), 
          alternative = "two.sided", adjust.method = c("BH", "holm", 
                                                       "hochberg", "hommel",
                                                       "bonferroni", "BY", "fdr", "none"), 
          plabel = c("p.adj.signif", "p.ad"), add = "jitter", xlab = NULL, 
          ylab = NULL, x.angle = 45, x.hjust = 1, x.vjust = 1, melt.data = T, 
          return.all = F, return.plot = F) {
  index = match.arg(arg = NULL, choices = index)
  plabel = match.arg(arg = NULL, choices = plabel)
  adjust.method = match.arg(arg = NULL, choices = adjust.method)
  df <- vegan_diver(meta, annotation = annotation, sample = sample, 
                    group = group, index = index, melt.data = melt.data)
  pairwise_t <- pairwise_t_test(df, adjust.method = adjust.method, 
                                alternative = alternative)
  if (is.null(ylab)) {
    ylab <- paste0(index, " index")
  }
  p <- plot_diver(df, pairwise_t, plabel = plabel, add = add, 
                  xlab = xlab, ylab = ylab)
  p <- p + theme(axis.text.x = element_text(angle = x.angle, 
                                            hjust = x.hjust, vjust = x.vjust))
  if (return.all) {
    return(list(df = df, pairwise_t = pairwise_t, p = p))
  }
  if (return.plot) {
    return(p)
  }
  return(pairwise_t)
}

#' T-test
#' 
#' @export
pairwise_t_test <- function (df, alternative = "two.sided", adjust.method = "BH") {
  pairwise_t <- df %>% rstatix::t_test(div_index ~ group, p.adjust.method = adjust.method, 
                                       alternative = alternative)
  return(pairwise_t)
}

#' boxplot
#' 
#' @export
plot_diver <- function (df, pairwise_t, plabel = c("p.adj.signif", "p.adj"), 
                        add = NULL, xlab = NULL, ylab = NULL) 
{
  plabel = match.arg(arg = NULL, choices = plabel)
  if (length(unique(df$group)) == 2) {
    plabel = NULL
  }
  p <- ggboxplot(df, x = "group", y = "div_index", add = add) + 
    stat_pvalue_manual(pairwise_t, label = plabel, y.position = seq(1:nrow(pairwise_t)) * 
                         max(df$div_index)/20 + max(df$div_index)) + xlab(xlab) + 
    ylab(ylab)
  p
}

#' Statistic of Group-enrichment
#'
#' @param x Cell clusters.
#' @param y Cell groups.
#' @param data data contain x and y. Default = NULL
#' @examples
#' data('meta200')
#' rovemat <- rove(x=as.vector(meta200$Annotation),y=as.vector(meta200$Group),plot=TRUE)
#' @export

rove <- function(data=NULL, x, y, plot=FALSE, normalize=FALSE,display_numbers = F,...){
  if(is.null(data)){
    df <- tibble(x=x, y=y)
  }else{
    df <- tibble(x=data[[x]], y=data[[y]])
  }
  n <- nrow(df)
  cont <- table(df)
  rs <- rowSums(cont)
  cs <- colSums(cont)
  basem <- n / as.matrix(rs) %*% t(as.matrix(cs)) * as.matrix.data.frame(cont)
  if(normalize)
    basem <- (basem - 1) * as.matrix(rep(1, length(rs))) %*% t(as.matrix(cs)) / n
  if(plot)
    pheatmap::pheatmap(basem, display_numbers = display_numbers,...)
  basem
}

#' Shannon, Simpson, and Fisher diversity indices and species richness.
#' 
#' @export
vegan_diver <- function (meta, annotation = "Annotation", sample = "Sample", 
          group = "Group", index = "simpson", melt.data = T) 
{
  allgroup <- as.vector(unique(meta[, group]))
  dflist <- lapply(allgroup, function(g) {
    mat <- table(meta[meta[, group] == g, c(sample, annotation)])
    simp <- vegan::diversity(mat, index = index)
    simdf <- data.frame(simp = simp, sample = names(simp), 
                        group = g, stringsAsFactors = F)
    colnames(simdf) <- c("div_index", "sample", "group")
    simdf
  })
  df <- do.call(rbind, dflist)
  if (!is.null(levels(meta[, group]))) {
    df$group <- factor(df$group, levels = levels(meta[, group]))
  }
  df
}

#' MASC - Mixed effect modeling of Associations of Single Cells
#'
#' @param dataset A data frame containing the contrast factor, random, and fixed effects for the model
#' @param cluster A factor indicating cluster assignments for each cell
#' @param contrast A vector indicating the variable to be tested for association with cluster abundance. Must match a column in dataset.
#' @param random_effects A vector indicating which terms should be modeled as random effects covariates. Terms listed must match columns in dataset.
#' @param fixed_effects A vector indicating which terms should be modeled as fixed effects covariates. Terms listed must match columns in dataset.
#' @param save_models Should MASC save the mixed-effects model objects generated for each cluster?
#' @param save_model_dir Location to save mixed-effect model objects. Defaults to current working directory.
#' @param verbose TRUE/FALSE
#'
#' @return data frame containing calculated association p-values and odds ratios for each cluster tested
#'
#' @author https://github.com/immunogenomics/masc
#' 
#' @examples
#' # Create test dataset with three clusters of 100 cells each
#' test.df <- data.frame(cluster = factor(rep(c(1, 2, 3), each = 100)))
#' # Create 6 donors that are cases or controls and include covariates
#' donors.df <- data.frame(donor = rep(paste("Donor", LETTERS[1:6], sep = "_"), each = 50),
#' sex = rep(c("M", "F", "M", "F", "F", "M"), each = 50),
#' status = rep(c("Case", "Case", "Control", "Control", "Case", "Control"), each = 50))
#' # Now make cluster 1 mostly case, cluster 2 mostly controls, etc
#' cases <- donors.df[donors.df$status == "Case",]
#' cases <- cases[sample(nrow(cases)),]
#' controls <- donors.df[donors.df$status == "Control",]
#' controls <- controls[sample(nrow(controls)),]
#' test.df <- cbind(rbind(cases[1:75,], controls[1:25,], cases[76:115,], controls[26:85,], cases[116:150,], controls[86:150,]), test.df)
#' # Test set call
#' library(lme4)
#' MASC(data = test.df, cluster = test.df$cluster, contrast = "status", random_effects = "donor", fixed_effects = "sex")
#'
#' @export
MASC <- function (dataset, cluster, contrast, random_effects = NULL, 
    fixed_effects = NULL, verbose = FALSE, save_models = FALSE, 
    save_model_dir = NULL) 
{
    if (is.factor(dataset[[contrast]]) == FALSE) {
        stop("Specified contrast term is not coded as a factor in dataset")
    }
    cluster <- as.character(cluster)
    designmat <- model.matrix(~cluster + 0, data.frame(cluster = cluster))
    dataset <- cbind(designmat, dataset)
    cluster <- as.character(cluster)
    designmat <- model.matrix(~cluster + 0, data.frame(cluster = cluster))
    dataset <- cbind(designmat, dataset)
    res <- vector(mode = "list", length = length(unique(cluster)))
    names(res) <- attributes(designmat)$dimnames[[2]]
    if (!is.null(fixed_effects) && !is.null(random_effects)) {
        model_rhs <- paste0(c(paste0(fixed_effects, collapse = " + "), 
            paste0("(1|", random_effects, ")", collapse = " + ")), 
            collapse = " + ")
        if (verbose == TRUE) {
            message(paste("Using null model:", "cluster ~", model_rhs))
        }
    }
    else if (!is.null(fixed_effects) && is.null(random_effects)) {
        model_rhs <- paste0(fixed_effects, collapse = " + ")
        if (verbose == TRUE) {
            message(paste("Using null model:", "cluster ~", model_rhs))
            stop("No random effects specified")
        }
    }
    else if (is.null(fixed_effects) && !is.null(random_effects)) {
        model_rhs <- paste0("(1|", random_effects, ")", collapse = " + ")
        if (verbose == TRUE) {
            message(paste("Using null model:", "cluster ~", model_rhs))
        }
    }
    else {
        model_rhs <- "1"
        if (verbose == TRUE) {
            message(paste("Using null model:", "cluster ~", model_rhs))
            stop("No random or fixed effects specified")
        }
    }
    cluster_models <- vector(mode = "list", 
                             length = length(attributes(designmat)$dimnames[[2]]))
    names(cluster_models) <- attributes(designmat)$dimnames[[2]]
    for (i in seq_along(attributes(designmat)$dimnames[[2]])) {
        test_cluster <- attributes(designmat)$dimnames[[2]][i]
        if (verbose == TRUE) {
            message(paste("Creating logistic mixed models for", 
                test_cluster))
        }
        null_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ 1 + "), 
            model_rhs), collapse = ""))
        full_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ ", 
            contrast, " + "), model_rhs), collapse = ""))
        null_model <- lme4::glmer(formula = null_fm, data = dataset, 
            family = binomial, nAGQ = 1, verbose = 0,
            control = lme4::glmerControl(optimizer = "bobyqa"))
        full_model <- lme4::glmer(formula = full_fm, data = dataset, 
            family = binomial, nAGQ = 1, verbose = 0, 
            control = lme4::glmerControl(optimizer = "bobyqa"))
        model_lrt <- anova(null_model, full_model)
        contrast_lvl2 <- paste0(contrast, levels(dataset[[contrast]])[2])
        contrast_ci <- lme4::confint.merMod(full_model, method = "Wald", 
            parm = contrast_lvl2)
        cluster_models[[i]]$null_model <- null_model
        cluster_models[[i]]$full_model <- full_model
        cluster_models[[i]]$model_lrt <- model_lrt
        cluster_models[[i]]$confint <- contrast_ci
    }
    output <- data.frame(cluster = attributes(designmat)$dimnames[[2]], 
        size = colSums(designmat))
    output$model.pvalue <- sapply(cluster_models, function(x) x$model_lrt[["Pr(>Chisq)"]][2])
    output[[paste(contrast_lvl2, "OR", sep = ".")]] <- sapply(cluster_models, 
        function(x) exp(lme4::fixef(x$full)[[contrast_lvl2]]))
    output[[paste(contrast_lvl2, "OR", "95pct.ci.lower", sep = ".")]] <- sapply(cluster_models, 
        function(x) exp(x$confint[contrast_lvl2, "2.5 %"]))
    output[[paste(contrast_lvl2, "OR", "95pct.ci.upper", sep = ".")]] <- sapply(cluster_models, 
        function(x) exp(x$confint[contrast_lvl2, "97.5 %"]))
    if (save_models == TRUE) {
        saveModelObj(cluster_models, save_dir = save_model_dir)
        return(output)
    }
    else {
        return(output)
    }
}

#' MASC - Mixed effect modeling of Associations of Single Cells
#' 
#' @examples
#' stat2x2 <- MASC_multi(meta200, cluster = 'maingroups', contrast = 'Group',
#' random_effects = 'Sample',model = '2x2')
#' stat_pairs <- MASC_multi(meta200, cluster = 'maingroups', contrast = 'Group',
#' random_effects = 'Sample',model = 'pairs')
#' @export
MASC_multi <- function (dataset, cluster, contrast, model = c("2x2", "pairs"), 
                        random_effects = NULL, fixed_effects = NULL, 
                        adjust.method = 'BH', verbose = FALSE, 
                        save_models = FALSE, save_model_dir = NULL, ncores = 4) 
{
  if (any(apply(dataset, 2, function(x) grepl(" ", x)))) {
    message("Warnning: Label with character ' ' is not allowed, and will be temply converted to '_space_'.")
    dataset <- as.data.frame(apply(dataset, 2, function(x) gsub(" ","_space_", x)))
  }
  if (any(apply(dataset, 2, function(x) grepl("-", x)))) {
    message("Warnning: Label with character '-' is not allowed, and will be temply converted to '_link_'.")
    dataset <- as.data.frame(apply(dataset, 2, function(x) gsub("-", "_link_", x)))
  }
  if (class(dataset[, contrast]) != "factor") {
    #stop(paste0("Column ", contrast, " in dataset must be factor."))
    dataset[, contrast] <- factor(dataset[, contrast], levels = sort(unique(dataset[, contrast])))
  }
  model <- match.arg(arg = NULL, choices = model)
  if (model == "2x2") {
    all_contrast <- levels(dataset[, contrast])
    dataset[, contrast] <- as.vector(dataset[, contrast])
    datalist <- parallel::mclapply(all_contrast, function(x) {
      dataset$status <- "other_status_as_control"
      dataset$status[dataset[, contrast] == x] <- x
      dataset$status <- factor(dataset$status,
                               levels = c("other_status_as_control", x))
      MASC_result <- MASC(dataset = dataset, cluster = dataset[, cluster], 
                          contrast = "status", random_effects = random_effects, 
                          fixed_effects = fixed_effects, verbose = verbose, 
                          save_models = save_models, save_model_dir = save_model_dir)
      result <- data.frame(Celltype = MASC_result[, "cluster"], 
                           Pvalue = MASC_result[, "model.pvalue"], 
                           OR = MASC_result[,  paste0("status", x, ".OR")], 
                           stringsAsFactors = F)
      result$Status <- x
      result$Celltype <- gsub("^cluster", "", result$Celltype)
      result$adjust.P <- p.adjust(result$Pvalue, method = adjust.method)
      return(result)
    }, mc.cores = ncores)
  }
  if (model == "pairs") {
    all_contrast <- levels(dataset[, contrast])
    dataset$status <- dataset[, contrast]
    dataset$status <- as.vector(dataset$status)
    all_pairs <- data.frame(i = rep(all_contrast, each = length(all_contrast)), 
                            j = rep(all_contrast, length(all_contrast)), 
                            stringsAsFactors = F)
    all_pairs <- all_pairs[all_pairs$i != all_pairs$j, ]
    datalist <- parallel::mclapply(1:nrow(all_pairs), function(x) {
      g1 <- as.character(all_pairs[x, "i"])
      g2 <- as.character(all_pairs[x, "j"])
      dataset <- dataset[dataset$status %in% c(g1, g2), ]
      if (length(unique(dataset$status)) < 2) next
      # g1 referred as control, g2 as query
      dataset$status <- factor(dataset$status, levels = c(g1, g2))
      MASC_result <- MASC(dataset = dataset, cluster = dataset[, cluster], 
                          contrast = "status", random_effects = random_effects, 
                          fixed_effects = fixed_effects, verbose = verbose, 
                          save_models = save_models, save_model_dir = save_model_dir)
      #print(head(MASC_result))
      result <- data.frame(Celltype = MASC_result[, "cluster"], 
                           Pvalue = MASC_result[, "model.pvalue"], 
                           OR = MASC_result[, paste0("status", g2, ".OR")],
                           stringsAsFactors = F)
      result$g1 <- g1
      result$g2 <- g2
      result$Celltype <- gsub("^cluster", "", result$Celltype)
      result$adjust.P <- p.adjust(result$Pvalue, method = adjust.method)
      return(result)
    }, mc.cores = ncores)
  }
  allresult <- do.call(rbind, datalist)
  allresult$Celltype <- gsub("_space_", " ", allresult$Celltype)
  allresult$Celltype <- gsub("_link_", "-", allresult$Celltype)
  return(allresult)
}


#' T-test of cell number
#' 
#' @export
cell_t_test <- function (metadf, Group = "Group", Annotation = "Annotation", 
          Sample = "Sample", alternative = c("two.sided", "great", 
                                             "less")) 
{
  alternative <- match.arg(arg = NULL, choices = alternative)
  metadf <- metadf[, c(Group, Annotation, Sample)]
  colnames(metadf) <- c("Group", "Annotation", "Sample")
  samdf <- metadf %>% group_by(Group, Annotation, Sample) %>% 
    summarise(count = n(), .groups = "drop")
  t_test_list <- lapply(unique(samdf$Annotation), function(cell) {
    dataset <- samdf[samdf$Annotation == cell, ]
    allgroup <- as.vector(unique(dataset$Group))
    t_test <- sapply(allgroup, function(x) {
      datag1 <- dataset[dataset$Group == x, "count", drop = T]
      datag2 <- dataset[dataset$Group != x, "count", drop = T]
      if (length(datag1) < 3 | length(datag2) < 3) {
        message("Warning: Not enough observations for ", 
                cell, " ", x, ". Return NA.")
        return(NA)
      }
      else {
        test <- t.test(datag1, datag2, alternative = alternative)
        test$p.value
      }
    })
    names(t_test) <- allgroup
    t_test
  })
  names(t_test_list) <- unique(samdf$Annotation)
  ttmat <- do.call(rbind, t_test_list)
  return(ttmat)
}


#' Statistics of cell proportion
#' 
#' @examples 
#' # 2x2
#' statdf <- cell_test(meta200, Group = 'Group',Annotation = 'Annotation',
#' Sample = 'Sample', model = '2x2', test = 'fisher',alternative = 'two.sided',
#' adjust.method = 'BH')
#' # pairwise
#' statdf <- cell_test(meta200, Group = 'Group',Annotation = 'Annotation',
#' Sample = 'Sample', model = '2x2', test = 'fisher',alternative = 'two.sided',
#' adjust.method = 'BH')
#' 
#' @export
cell_test <- function (dataset, Group = "Group", Annotation = "Annotation", 
                       Sample = "Sample", 
                       test = c("fisher", "chisq", "t", "wilcox"), 
                       model = c("2x2", "pairs"), 
                       alternative = c("two.sided", "great", "less"), 
                       adjust.method = c("holm", "hochberg", "hommel",
                                         "bonferroni", "BH", "BY", "fdr",
                                         "none")) {
    test <- match.arg(arg = NULL, choices = test)
    alternative <- match.arg(arg = NULL, choices = alternative)
    metadf <- data.frame(Group = as.vector(dataset[, Group]), 
        Annotation = as.vector(dataset[, Annotation]), Sample = as.vector(dataset[, 
            Sample]), stringsAsFactors = F)
    samdf <- metadf %>% group_by(Group, Annotation, Sample) %>% 
        summarise(count = n(), .groups = "drop")
    if (model == "2x2") {
        fisher_test_list <- lapply(unique(samdf$Annotation), function(cell) {
                datacell <- samdf
                datacell[datacell$Annotation != cell, "Annotation"] <- "other_cell_as_control"
                allgroup <- as.vector(unique(samdf$Group))
                fisher_test_result <- lapply(allgroup, function(x) {
                  datag <- datacell
                  datag[datag$Group != x, "Group"] <- "other_group_as_control"
                  summat <- datag %>% group_by(Group, Annotation) %>% 
                    summarise(sum = sum(count), .groups = "drop")
                  matd <- reshape2::dcast(summat, Group ~ Annotation, 
                    value.var = "sum")
                  matd <- data.frame(matd, row.names = 1, check.names = F)
                  matd <- matd[c(x, "other_group_as_control"), 
                    c(cell, "other_cell_as_control")]
                  matd[is.na(matd)] <- 0
                  ft <- rstatix::fisher_test(matd, detailed = TRUE,
                                             alternative = alternative)
                  OR_p <- c(ft$estimate, ft$p)
                  if (test == "chisq") {
                    chisq <- chisq.test(matd)
                    OR_p[2] <- chisq$p.value
                  }
                  names(OR_p) <- c("OR", "Pvalue")
                  return(OR_p)
                })
                names(fisher_test_result) <- allgroup
                ftmat <- as.data.frame(do.call(rbind, fisher_test_result))
                ftmat$Annotation <- cell
                ftmat$Group <- rownames(ftmat)
                return(ftmat)
            })
        names(fisher_test_list) <- unique(samdf$Annotation)
        ftall <- do.call(rbind, fisher_test_list)
        ftall$adj.Pvalue <- p.adjust(ftall$Pvalue, method = adjust.method)
    }
    if (model == "pairs") {
        fisher_test_list <- lapply(unique(samdf$Annotation), 
            function(cell) {
                datacell <- samdf
                datacell[datacell$Annotation != cell, "Annotation"] <- "other_cell_as_control"
                summat <- datacell %>% group_by(Group, Annotation) %>% 
                  summarise(sum = sum(count), .groups = "drop")
                datag <- reshape2::dcast(summat, formula = Annotation ~ 
                  Group, value.var = "sum")
                datag[is.na(datag)] <- 0
                datag <- data.frame(datag, row.names = 1)
                datag <- datag[c(cell, "other_cell_as_control"), 
                  ]
                OR_p <- rstatix::pairwise_fisher_test(datag, 
                  detailed = TRUE, alternative = alternative, 
                  p.adjust.method = adjust.method)
                if (test == "chisq") {
                  OR_p <- rstatix::pairwise_chisq_gof_test(datag, 
                    p.adjust.method = adjust.method)
                }
                OR_p$Annotation <- cell
                return(OR_p)
            })
        names(fisher_test_list) <- unique(samdf$Annotation)
        ftall <- do.call(rbind, fisher_test_list)
    }
    return(ftall)
}

