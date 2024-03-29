# the main script is from https://github.com/waldronlab/lefser/blob/master/R/lefser.R
fillPmatZmat <- function(group,
                         block,
                         expr_sub,
                         p.threshold){
  # creates a list of boolean vectors, each vector indicates
  # existance (TRUE) or absence (FALSE) of a class/sub-class combination
  combos <- apply(
    expand.grid(levels(group), levels(block)), 1L, paste0, collapse = "")
  combined <- paste0(as.character(group), as.character(block))
  logilist <- lapply(setNames(nm = sort(combos)), `==`, combined)

  ## uses Wilcoxon rank-sum test to test for significant differential abundances between
  ## subclasses of one class against subclasses of all othe classes; results are saved in
  ## "pval_mat" and "z_mat" matrices
  whichlist <- lapply(logilist, which)
  sblock <- seq_along(levels(block))
  texp_sub <- t(expr_sub)
  iters <- expand.grid(sblock, sblock + length(sblock))
  group_formats <- apply(iters, 1L, function(x) {
    ind <- unlist(whichlist[x])
    apply(texp_sub, 2L, function(g) {
      wx <- suppressWarnings(coin::wilcox_test(g ~ group, subset = ind))
      cbind.data.frame(
        p.value = coin::pvalue(wx), statistic = coin::statistic(wx)
      )
    })
  })

  res <- lapply(group_formats, function(x) do.call(rbind, x))
  pval_mat <- do.call(cbind, lapply(res, `[[`, "p.value"))
  z_mat <- do.call(cbind, lapply(res, `[[`, "statistic"))

  rownames(pval_mat) <- rownames(expr_sub)
  rownames(z_mat) <- rownames(expr_sub)

  ## converts "pval_mat" into boolean matrix "logical_pval_mat" where
  ## p-values <= wilcoxon.threshold
  logical_pval_mat <- pval_mat <= p.threshold * 2.0
  logical_pval_mat[is.na(logical_pval_mat)] <- FALSE

  ## determines which rows (features) have all p-values<=0.05
  ## and selects such rows from the matrix of z-statistics
  sub <- apply(logical_pval_mat, 1L, all)
  z_mat_sub <- z_mat[sub, , drop = FALSE]
  # confirms that z-statistics of a row all have the same sign
  sub <- abs(rowSums(z_mat_sub)) == rowSums(abs(z_mat_sub))
  expr_sub[names(sub[sub]), , drop = FALSE]
}

## ensures that more than half of the values in each for each feature are unique
## if that is not the case then a count value is altered by adding it to a small value
## generated via normal distribution with mean=0 and sd=5% of the count value
createUniqueValues <- function(df, group){
  orderedrows <- rownames(df)
  splitdf <- split(df, group)
  maxim <- vapply(table(group), function(x) max(x * 0.5, 4), numeric(1L))
  for (i in seq_along(splitdf)) {
    sdat <- splitdf[[i]]
    splitdf[[i]][] <- lapply(sdat, function(cols) {
      if (length(unique(cols)) > maxim[i])
        cols
      else
        abs(cols + stats::rnorm(
          length(cols), mean = 0, sd = max(cols * 0.05, 0.01))
        )
    })
  }
  df <- do.call(rbind, unname(splitdf))
  df[match(orderedrows, rownames(df)),, drop = FALSE]
}

contastWithinClassesOrFewPerClass <- function(expr_sub_t_df, rand_s, min_cl, ncl, groups){

  cols <- expr_sub_t_df[rand_s, , drop = FALSE]
  cls <- expr_sub_t_df$class[rand_s]
  # if the number of classes is less than the actual number (typically two)
  # of classes in the dataframe then return TRUE
  if (length(unique(cls)) < ncl) {
    return (TRUE)
  }
  # detect if for each class there are not fewer than the minimum (min_cl) number of samples
  if (any(table(cls) < min_cl)) {
    return (TRUE)
  }
  # separate the randomly selected samples (cols) into a list of the two classes
  drops <- c("class")
  by_class <-
      lapply(seq_along(groups), function(x) {
        cols[cols[, "class"] == groups[x],!(names(cols) %in% drops)]
    })

  # makes sure that within each class all features have at least min_cl unique count values
  for (i in seq_along(groups)) {
      unique_counts_per_microb = apply(by_class[[i]], 2, function(x) {
        length(unique(x))
      })
      if ((any(unique_counts_per_microb <= min_cl) &
           min_cl > 1) |
          (min_cl == 1 & any(unique_counts_per_microb <= 1))) {
        return (TRUE)
      }
    }
    return (FALSE)

}

ldaFunction <- function (data, lfk, rfk, min_cl, ncl, groups){
  # test 1000 samples for contrast within classes per feature
  # and that there is at least a minimum number of samples per class
  for (j in 1:1000) {
    rand_s <- sample(seq_len(lfk), rfk, replace = TRUE)
    if (!contastWithinClassesOrFewPerClass(data, rand_s, min_cl, ncl, groups)) {
      break
    }
  }
  # lda with rfk number of samples
  lda.fit <- MASS::lda(class ~ ., data = data, subset = rand_s)
  # coefficients that transform observations to discriminants
  w <- lda.fit$scaling[, 1]
  # scaling of lda coefficients
  w.unit <- w / sqrt(sum(w ^ 2))
  sub_d <- data[rand_s,]
  ss <- sub_d[,-match("class", colnames(sub_d))]
  xy.matrix <- as.matrix(ss)
  # the original matrix is transformed
  LD <- xy.matrix %*% w.unit
  # effect size is calculated as difference between averaged disciminants
  # of two classes
  effect_size <-
    abs(mean(LD[sub_d[, "class"] == 0]) - mean(LD[sub_d[, "class"] == 1]))
  # scaling lda coefficients by the efect size
  scal <- w.unit * effect_size
  # mean count values per fclass per feature
  rres <- lda.fit$means
  rowns <- rownames(rres)
  lenc <- length(colnames(rres))

  coeff <- vector("numeric", length(scal))
  for (v in seq_along(scal)) {
    if (!is.na(scal[v])) {
      coeff[v] <- abs(scal[v])
    } else{
      coeff[v] <- 0
    }
  }
  # count value differences between means of two classes for each feature
  lda.means.diff <- (lda.fit$means[2,] - lda.fit$means[1,])
  # difference between a feature's class means and effect size adjusted lda coefficient
  # are averaged for each feature
  (lda.means.diff + coeff) / 2
}

.numeric01 <- function(x) {
  x <- as.factor(x)
  uvals <- levels(x)
  ifelse(x == uvals[1L], 0L, 1L)
}

filterKruskal <- function(expr, group, p.value) {
  # applies "kruskal.test.alt" function to each row (feature) of expr
  # to detect differential abundance between classes, 0 and 1
  kw.res <- apply(expr, 1L, function(x) {
    stats::kruskal.test(x ~ group)[["p.value"]]
  })
  # selects p-values less than or equal to kw.threshold
  kw.sub <- kw.res <= p.value

  # eliminates NAs
  kw.sub[is.na(kw.sub)] <- FALSE

  # extracts features with statistically significant differential abundance
  # from "expr" matrix
  expr[kw.sub,]
}

.trunc <- function(scores_df, trim.names){
  Names <- gsub("`", "", scores_df[["Names"]])
  if (trim.names) {
    listNames <- strsplit(Names, "\\||\\.")
    Names <- vapply(listNames, tail, character(1L), 1L)
  }
  scores_df[["Names"]] <- Names
  return(scores_df)
}

#' @title R implementation of the LEfSe method
#'
#' @description
#' Perform a LEfSe analysis: the function carries out differential analysis
#' between two sample groups for multiple microorganisms and uses linear discirminant analysis
#' to establish their effect sizes. Subclass information for each class can be incorporated
#' into the analysis (see examples). Microorganisms with large differences between two sample groups
#' are identified as biomarkers.
#'
#' @param expr A \code{\linkS4class{SummarizedExperiment}} with expression data.
#' @param kruskal.threshold numeric(1) The p-value for the Kruskal-Wallis Rank
#' Sum Test (default 0.05).
#' @param wilcox.threshold numeric(1) The p-value for the Wilcoxon Rank-Sum Test
#' when 'blockCol' is present (default 0.05).
#' @param lda.threshold numeric(1) The effect size threshold (default 2.0).
#' @param groupCol character(1) Column name in `colData(expr)` indicating
#' groups, usually a factor with two levels (e.g., `c("cases", "controls")`;
#' default "GROUP").
#' @param blockCol character(1) Optional column name in `colData(expr)`
#' indicating the blocks, usually a factor with two levels (e.g.,
#' `c("adult", "senior")`; default NULL).
#' @param assay The i-th assay matrix in the `SummarizedExperiment` ('expr';
#' default 1).
#' @param trim.names If `TRUE` extracts the most specific taxonomic rank of organism.
#' @return
#' The function returns a dataframe with two columns, which are
#' names of microorganisms and their LDA scores.
#'
#' @export
#'
#' @importFrom stats kruskal.test reorder rnorm
#' @importFrom coin pvalue statistic wilcox_test
#' @importFrom MASS lda
#' @importFrom methods as is
#' @importFrom stats setNames
#' @importFrom utils tail
#' @import SummarizedExperiment
#'
#' @examples
#'
lefser <- function(expr,
           kruskal.threshold = 0.05,
           wilcox.threshold = 0.05,
           lda.threshold = 2.0,
           groupCol = "GROUP",
           blockCol = NULL,
           assay = 1L,
           trim.names = FALSE)
  {
    groupf <- colData(expr)[[groupCol]]
    if (is.null(groupf))
      stop("A valid group assignment 'groupCol' must be provided")
    groupf <- as.factor(groupf)
    groupsf <- levels(groupf)
    if (length(groupsf) != 2L)
      stop(
        "Group classification is not dichotomous:\n",
        "Found (", paste(groupsf, collapse = ", "), ")"
      )
    group <- .numeric01(groupf)
    groups <- 0:1
    expr_data <- assay(expr, i = assay)
    expr_sub <- filterKruskal(expr_data, group, kruskal.threshold)

    if (!is.null(blockCol)) {
      block <- as.factor(colData(expr)[[blockCol]])
      expr_sub <- fillPmatZmat(groupf, block, expr_sub, wilcox.threshold)
    }

    # transposes matrix and add a "class" (i.e., group) column
    # matrix converted to dataframe
    expr_sub_t <- t(expr_sub)
    expr_sub_t_df <- as.data.frame(expr_sub_t)
    expr_sub_t_df <- createUniqueValues(expr_sub_t_df, groupf)
    expr_sub_t_df <- cbind(expr_sub_t_df, class = group)

    # number of samples (i.e., subjects) in the dataframe
    lfk <- nrow(expr_sub_t_df)
    # rfk is the number of subject that will be used in linear discriminant analysis
    rfk <- floor(lfk * 2 / 3)
    # number of classes (two)
    ncl <- length(groups)
    # count samples in each class of the dataframe, select the number from the class with a smaller
    # count of samples and multiply that number by 2/*2/3*0.5
    min_cl <-
      as.integer(min(table(expr_sub_t_df$class)) * 2 / 3 * 2 / 3 *
                   0.5)
    # if min_cl is less than 1, then make it equal to 1
    min_cl <- max(min_cl, 1)

    # lda_fn repeated 30 times, producing a matrix of 30 scores per feature
    eff_size_mat <-
      replicate(30, suppressWarnings(ldaFunction(
        expr_sub_t_df, lfk, rfk, min_cl, ncl, groups
      )), simplify = TRUE)

    # mean of 30 scores per feature
    raw_lda_scores <- rowMeans(eff_size_mat)

    # processing of score
    processed_scores <-
      sign(raw_lda_scores) * log((1 + abs(raw_lda_scores)), 10)

    # sorting of scores
    processed_sorted_scores <- sort(processed_scores)
    scores_df <- data.frame(Names = names(processed_sorted_scores),
                            scores = as.vector(processed_sorted_scores),
                            stringsAsFactors = FALSE)

    scores_df <- .trunc(scores_df, trim.names)

    threshold_scores <- abs(scores_df$scores) >= lda.threshold
    scores_df[threshold_scores, ]
  }

############################################################################
#' @title Differential Expression Analysis by LEfSe Package
#'
#' @description
#' LEfSe method for microbiome biomarker discovery.
#'
#' @details 12/2/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param Expression, ExpressionSet; (Required) ExpressionSet object.
#' @param trim, Character; filter to apply.(default: trim="none").
#' @param transform, Character; transformation to apply.(default: tranform="none").
#' @param normalize, Character; normalization to apply.(default: normalize="none").
#' @param Group_info, Character; design factor(default: "Group").
#' @param Group_name, Character; (Required) the group for comparison.
#' @param kw.p, Numeric; significant level of kruskal.test(default: 0.05).
#' @param wl.p, Numeric; significant level of wilcox(default: 0.05).
#' @param Lda, Numeric; LDA score(default: 2).
#'
#' @return
#'   significant difference with enriched directors
#'
#' @export
#'
#' @importFrom dplyr %>% select filter intersect all_of mutate
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom stats setNames sd
#' @import SummarizedExperiment
#' @importFrom Biobase pData exprs
#' @importFrom MASS lda
#'
#' @usage run_LEfSe(dataset=ExpressionSet,
#'                  trim="none",
#'                  transform="none",
#'                  normalize="none",
#'                  Group_info="Group",
#'                  Group_name=c("HC", "AA"),
#'                  kw.p=0.05,
#'                  wl.p=0.05,
#'                  Lda=2)
#' @examples
#'
#' \donttest{
#' data(ExprSetRawCount)
#'
#' LEfSe_res <- run_LEfSe(dataset=ExprSetRawCount, Group_info="Group", Group_name=c("HC", "AA"), kw.p=0.05, wl.p=0.05, Lda=2)
#' LEfSe_res
#' }
#'
run_LEfSe <- function(dataset=ExprSetRawCount,
                      trim="none",
                      transform="none",
                      normalize="none",
                      Group_info="Group",
                      Group_name=c("HC", "AA"),
                      kw.p=0.05,
                      wl.p=0.05,
                      Lda=2){

  # preprocess
  dataset_processed <- get_processedExprSet(dataset=dataset,
                                            trim=trim,
                                            transform=transform,
                                            normalize=normalize)

  metadata <- Biobase::pData(dataset_processed)
  colnames(metadata)[which(colnames(metadata) == Group_info)] <- "Group"
  profile <- Biobase::exprs(dataset_processed)
  if(!any(profile %% 1 == 0)){
    stop("The input matrix is not integer matrix please Check it")
  }

  # choose group
  phen <- metadata %>% dplyr::filter(Group%in%Group_name)
  intersect_sid <- dplyr::intersect(rownames(phen), colnames(profile))
  # Prepare for input data
  colData <- phen %>% tibble::rownames_to_column("SampleID") %>%
    dplyr::select(dplyr::all_of(c("SampleID", "Group"))) %>%
    dplyr::filter(SampleID%in%intersect_sid) %>%
    dplyr::mutate(Group=factor(Group, levels = Group_name)) %>%
    tibble::column_to_rownames("SampleID")
  proData <- profile %>% data.frame() %>%
    dplyr::select(dplyr::all_of(rownames(colData))) %>%
    as.matrix()

  if(!all(rownames(colData) == colnames(proData))){
    stop("Order of sampleID between colData and proData is wrong please check your data")
  }

  # Median abundance
  median_res <- apply(proData, 1, function(x, y){
    dat <- data.frame(value=as.numeric(x), group=y)
    # median value
    mn <- tapply(dat$value, dat$group, median) %>%
      data.frame() %>% setNames("value") %>%
      tibble::rownames_to_column("Group")
    mn1 <- with(mn, mn[Group%in%Group_name[1], "value"])
    mn2 <- with(mn, mn[Group%in%Group_name[2], "value"])
    mnall <- median(dat$value)

    res <- c(mnall, mn1, mn2)
    return(res)
  }, colData$Group) %>%
    t() %>% data.frame() %>%
    tibble::rownames_to_column("FeatureID")
  colnames(median_res) <- c("FeatureID", "Median Abundance\n(All)",
                            paste0("Median Abundance\n", Group_name))

  se <- SummarizedExperiment::SummarizedExperiment(
              assays=list(counts=proData),
              colData=colData,
              metadata="Profile")
  lefse_res <- lefser(se,
                      kruskal.threshold = kw.p,
                      wilcox.threshold  = wl.p,
                      lda.threshold     = Lda,
                      groupCol = Group_info,
                      blockCol = NULL,
                      assay    = 1L,
                      trim.names = TRUE)

  t_res <- lefse_res %>% dplyr::mutate(Group=ifelse(scores > 0, Group_name[2], Group_name[1]))
  colnames(t_res) <- c("FeatureID", "LDA_Score", "Enrichment")

  t_res_temp <- t_res %>%
    dplyr::inner_join(median_res, by = "FeatureID")

  # Number of Group
  dat_status <- table(colData$Group)
  dat_status_number <- as.numeric(dat_status)
  dat_status_name <- names(dat_status)
  t_res_temp$Block <- paste(paste(dat_status_number[1], dat_status_name[1], sep = "_"),
                       "vs",
                       paste(dat_status_number[2], dat_status_name[2], sep = "_"))

  t_res_temp2 <- t_res_temp %>% dplyr::select(FeatureID, Block, Enrichment, LDA_Score,
                                              dplyr::everything()) %>% dplyr::arrange(LDA_Score)

  # 95% CI Odd Ratio
  res_odd <- run_OddRatio(colData, proData, Group_name)

  # Merge
  res <- dplyr::inner_join(t_res_temp2, res_odd, by="FeatureID")

  return(res)
}
