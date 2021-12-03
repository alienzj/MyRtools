# https://github.com/FrederickHuangLin/ANCOM/blob/master/scripts/ancom_v2.1.R
#' @title Differential Expression Analysis by Analysis of composition of microbiomes(ANCOM)
#'
#' @description
#' ANCOM accounts for the underlying structure in the data and can be used for comparing
#' the composition of microbiomes in two or more populations.
#' ANCOM makes no distributional assumptions and can be implemented in
#' a linear model framework to adjust for covariates as well as model longitudinal data.
#' ANCOM also scales well to compare samples involving thousands of taxa.
#'
#' @references Mandal et al. "Analysis of composition of microbiomes: a novel
#' method for studying microbial composition", Microbial Ecology in Health
#' & Disease, (2015), 26.
#'
#' @details 12/3/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param Expression, ExpressionSet; (Required) ExpressionSet object.
#' @param GroupVar, Character; design factor(default: "Group").
#' @param AdjVar, Character; the adjusted variables.
#' @param RandVar, Character; random effects for longitudinal analysis or repeated measure("~ 1 | studyid").
#' @param Pvalue, Numeric; significant level(default: 0.05).
#' @param Wvalue, Numeric; W statistic for clarify the significant features(default: 0.7).
#'
#' @return
#' a list of results:
#'   significant features pass the threshold of W statistics
#'   A volcano plot of significant features
#'
#' @importFrom dplyr %>% select filter intersect pull all_of
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom stats setNames lm residuals quantile dnorm sd wilcox.test kruskal.test aov na.omit formula
#' @importFrom compositions clr
#' @importFrom nlme lme
#' @importFrom tidyr gather
#' @import convert
#' @import ggplot2
#'
#' @usage DA_ANCOM(dataset=ExpressionSet, GroupVar="Group", AdjVar=NULL, RandVar=NULL, Pvalue=0.05, Wvalue=0.7)
#' @examples
#'
#' data(ExprSet_species_count)
#'
#' ANCOM_res <- DA_ANCOM(dataset=ExprSet_species_count, GroupVar="Group", AdjVar="Gender", RandVar=NULL, Pvalue=0.05, Wvalue=0.7)
#' ANCOM_res$res
#'
DA_ANCOM <- function(dataset=ExprSet_species_count,
                     GroupVar="Group",
                     AdjVar=NULL,
                     RandVar=NULL,
                     Pvalue=0.05,
                     Wvalue=0.7){

  # Data Pre-Processing
  preprocess_ANCOM <- function(
    feature_table,
    meta_data,
    sample_var = NULL,
    group_var = NULL,
    out_cut = 0.05,
    zero_cut = 0.90,
    lib_cut = 100,
    neg_lb = FALSE){

    feature_table <- data.frame(feature_table, check.names = FALSE)
    meta_data <- data.frame(meta_data, check.names = FALSE)
    # Drop unused levels
    meta_data[] <- lapply(
      meta_data,
      function(x)if(is.factor(x)) factor(x) else x
    )
    # Match sample IDs between metadata and feature table
    sample_ID <- dplyr::intersect(meta_data[, sample_var], colnames(feature_table))
    feature_table <- feature_table[, sample_ID]
    meta_data <- meta_data[match(sample_ID, meta_data[, sample_var]), ]

    # 1. Identify outliers within each taxon
    if(!is.null(group_var)){
      group <- meta_data[, group_var]
      z <- feature_table + 1 # Add pseudo-count (1)
      f <- log(z)
      f[f == 0] <- NA
      f <- colMeans(f, na.rm = TRUE)
      f_fit <- stats::lm(f ~ group)
      e <- rep(0, length(f))
      e[!is.na(group)] <- stats::residuals(f_fit)
      y <- t(t(z) - e)

      outlier_check <- function(x){
        # Fitting the mixture model using the algorithm of Peddada, S. Das, and JT Gene Hwang (2002)
        mu1 <- stats::quantile(x, 0.25, na.rm = TRUE)
        mu2 <- stats::quantile(x, 0.75, na.rm = TRUE)
        sigma1 <- stats::quantile(x, 0.75, na.rm = TRUE) - stats::quantile(x, 0.25, na.rm = TRUE)
        sigma2 <- sigma1
        pi <- 0.75
        n <- length(x)
        epsilon <- 100
        tol <- 1e-5
        score <- pi * stats::dnorm(x, mean = mu1, sd = sigma1)/
          ((1 - pi) * stats::dnorm(x, mean = mu2, sd = sigma2))
        while(epsilon > tol){
          grp1_ind <- (score >= 1)
          mu1_new <-  mean(x[grp1_ind])
          mu2_new <- mean(x[!grp1_ind])
          sigma1_new <- sd(x[grp1_ind])
          if(is.na(sigma1_new)){sigma1_new <- 0}
          sigma2_new <- stats::sd(x[!grp1_ind])
          if(is.na(sigma2_new)) {sigma2_new <- 0}
          pi_new <- sum(grp1_ind)/n

          para <- c(mu1_new, mu2_new, sigma1_new, sigma2_new, pi_new)
          if(any(is.na(para))) break

          score <- pi_new * stats::dnorm(x, mean = mu1_new, sd = sigma1_new)/
            ((1-pi_new) * stats::dnorm(x, mean = mu2_new, sd = sigma2_new))

          epsilon <- sqrt(
            (mu1 - mu1_new)^2 +
              (mu2 - mu2_new)^2 +
              (sigma1 - sigma1_new)^2 +
              (sigma2 - sigma2_new)^2 +
              (pi - pi_new)^2)
          mu1 <- mu1_new
          mu2 <- mu2_new
          sigma1 <- sigma1_new
          sigma2 <- sigma2_new
          pi <- pi_new
        }

        if(mu1 + 1.96 * sigma1 < mu2 - 1.96 * sigma2){
          if(pi < out_cut){
            out_ind <- grp1_ind
          }else if(pi > 1 - out_cut){
            out_ind <- (!grp1_ind)
          }else{
            out_ind <- rep(FALSE, n)
          }
        }else{
          out_ind <- rep(FALSE, n)
        }
        return(out_ind)
      }
      out_ind <- matrix(FALSE, nrow = nrow(feature_table), ncol = ncol(feature_table))
      out_ind[, !is.na(group)] <- t(apply(y, 1, function(i){
        unlist(tapply(i, group, function(j){outlier_check(j)}))
      }
      )
      )

      feature_table[out_ind] <- NA
    }

    # 2. Discard taxa with zeros  >=  zero_cut
    zero_prop <- apply(feature_table, 1, function(x){
      sum(x == 0, na.rm = TRUE)/length(x[!is.na(x)])
    }
    )
    taxa_del <- which(zero_prop >= zero_cut)
    if(length(taxa_del) > 0){
      feature_table <- feature_table[- taxa_del, ]
    }

    # 3. Discard samples with library size < lib_cut
    lib_size <- colSums(feature_table, na.rm = TRUE)
    if(any(lib_size < lib_cut)){
      subj_del <- which(lib_size < lib_cut)
      feature_table <- feature_table[, - subj_del]
      meta_data <- meta_data[- subj_del, ]
    }

    # 4. Identify taxa with structure zeros
    if(!is.null(group_var)){
      group <- factor(meta_data[, group_var])
      present_table <- as.matrix(feature_table)
      present_table[is.na(present_table)] <- 0
      present_table[present_table != 0] <- 1

      p_hat <- t(apply(present_table, 1, function(x){
        unlist(tapply(x, group, function(y) mean(y, na.rm = TRUE)))
          }
        )
      )
      samp_size <- t(apply(feature_table, 1, function(x){
        unlist(tapply(x, group, function(y) length(y[!is.na(y)])))
          }
        )
      )
      p_hat_lo <- p_hat - 1.96 * sqrt(p_hat * (1 - p_hat)/samp_size)

      struc_zero <- (p_hat == 0) * 1
      # Whether we need to classify a taxon into structural zero by its negative lower bound?
      if(neg_lb){
        struc_zero[p_hat_lo <= 0] <- 1
      }

      # Entries considered to be structural zeros are set to be 0s
      struc_ind <- struc_zero[, group]
      feature_table <- feature_table * (1 - struc_ind)

      colnames(struc_zero) <- paste0("structural_zero (", colnames(struc_zero), ")")
    }else{
      struc_zero <- NULL
    }

    # 5. Return results
    res <- list(feature_table=feature_table,
                meta_data=meta_data,
                structure_zeros=struc_zero)
    return(res)
  }

  # ANCOM main function
  ANCOM <- function(feature_table,
                    meta_data,
                    struc_zero = NULL,
                    group_var = NULL,
                    p_adj_method = "BH",
                    alpha = 0.05,
                    adj_formula = NULL,
                    rand_formula = NULL,
                    w_cutoff = Wvalue,
                    ...){

    # OTU table transformation:
    # (1) Discard taxa with structural zeros (if any); (2) Add pseudocount (1) and take logarithm.
    if (!is.null(struc_zero)) {
      num_struc_zero <- apply(struc_zero, 1, sum)
      comp_table <- feature_table[num_struc_zero == 0, ]
    }else{
      comp_table <- feature_table
    }
    comp_table <- log(as.matrix(comp_table) + 1)
    n_taxa <- dim(comp_table)[1]
    taxa_id <- rownames(comp_table)
    n_samp <- dim(comp_table)[2]

    # Determine the type of statistical test and its formula.
    if(is.null(rand_formula) & is.null(adj_formula)){
      # Basic model
      # Whether the main variable of interest has two levels or more?
      if(length(unique(meta_data %>% dplyr::pull(group_var))) == 2){
        # Two levels: Wilcoxon rank-sum test
        tfun <- stats::wilcox.test
      }else{
        # More than two levels: Kruskal-Wallis test
        tfun <- stats::kruskal.test
      }
      # Formula
      tformula <- stats::formula(paste("x ~", group_var, sep = " "))
    }else if(is.null(rand_formula) & !is.null(adj_formula)){
      # Model: ANOVA
      tfun <- stats::aov
      # Formula
      tformula <- stats::formula(paste("x ~", group_var, "+", paste(adj_formula, collapse = "+"), sep = " "))
    }else if(!is.null(rand_formula)){
      # Model: Mixed-effects model
      tfun <- nlme::lme
      # Formula
      if(is.null(adj_formula)){
        # Random intercept model
        tformula <- stats::formula(paste("x ~", group_var))
      }else{
        # Random coefficients/slope model
        tformula <- stats::formula(paste("x ~", group_var, "+", paste(adj_formula, collapse = "+")))
      }
    }

    # Calculate the p-value for each pairwise comparison of taxa.
    p_data <- matrix(NA, nrow = n_taxa, ncol = n_taxa)
    colnames(p_data) <- taxa_id
    rownames(p_data) <- taxa_id
    for(i in 1:(n_taxa - 1)){
      # Loop through each taxon.
      # For each taxon i, additive log ratio (alr) transform the OTU table using taxon i as the reference.
      # e.g. the first alr matrix will be the log abundance data (comp_table) recursively subtracted
      # by the log abundance of 1st taxon (1st column) column-wisely, and remove the first i columns since:
      # the first (i - 1) columns were calculated by previous iterations, and
      # the i^th column contains all zeros.
      alr_data <- apply(comp_table, 1, function(x) x - comp_table[i, ])
      # apply(...) allows crossing the data in a number of ways and avoid explicit use of loop constructs.
      # Here, we basically want to iteratively subtract each column of the comp_table by its i^th column.
      alr_data <- alr_data[, - (1:i), drop = FALSE]
      n_lr <- dim(alr_data)[2] # number of log-ratios (lr)
      alr_data <- cbind(alr_data, meta_data) # merge with the metadata

      # P-values
      if(is.null(rand_formula) & is.null(adj_formula)){
        p_data[-(1:i), i] <- apply(alr_data[, 1:n_lr, drop = FALSE], 2, function(x){
          suppressWarnings(tfun(tformula,
                                data = data.frame(x, alr_data,
                                                  check.names = FALSE))$p.value)
          }
        )
      }else if(is.null(rand_formula) & !is.null(adj_formula)) {
        p_data[-(1:i), i] <- apply(alr_data[, 1:n_lr, drop = FALSE], 2, function(x){
          fit <- tfun(tformula,
                      data = data.frame(x, alr_data, check.names = FALSE),
                      na.action = stats::na.omit)
          summary(fit)[[1]][group_var, "Pr(>F)"]
          }
        )
      }else if(!is.null(rand_formula)){
        p_data[-(1:i), i] <- apply(alr_data[, 1:n_lr, drop = FALSE], 2, function(x){
          fit <- tfun(fixed = tformula,
                      data = data.frame(x, alr_data, check.names = FALSE),
                      random = formula(rand_formula),
                      na.action = stats::na.omit, ...)
          anova(fit)[group_var, "p-value"]
          }
        )
      }
    }
    # Complete the p-value matrix.
    # What we got from above iterations is a lower triangle matrix of p-values.
    p_data[upper.tri(p_data)] <- t(p_data)[upper.tri(p_data)]
    diag(p_data) <- 1 # let p-values on diagonal equal to 1
    p_data[is.na(p_data)] <- 1 # let p-values of NA equal to 1

    # Multiple comparisons correction.
    q_data <- apply(p_data, 2, function(x) p.adjust(x, method = p_adj_method))

    # Calculate the W statistic of ANCOM.
    # For each taxon, count the number of q-values < alpha.
    W <- apply(q_data, 2, function(x) sum(x < alpha))

    # Organize outputs
    out_comp <- data.frame(taxa_id, W, row.names = NULL, check.names = FALSE)
    # Declare a taxon to be differentially abundant based on the quantile of W statistic.
    # We perform (n_taxa - 1) hypothesis testings on each taxon, so the maximum number of rejections is (n_taxa - 1).
    out_comp <- out_comp %>% dplyr::mutate(
      detected_Wcutoff = ifelse(W > w_cutoff * (n_taxa -1), TRUE, FALSE))

    # Taxa with structural zeros are automatically declared to be differentially abundant
    if(!is.null(struc_zero)){
      out <- data.frame(taxa_id = rownames(struc_zero),
                        W = Inf,
                        detected_Wcutoff = TRUE,
                        row.names = NULL,
                        check.names = FALSE)
      out[match(taxa_id, out$taxa_id), ] <- out_comp
    }else{
      out <- out_comp
    }
    colnames(out)[1:3] <- c("FeatureID", "W", paste0("detected_", w_cutoff))

    # Draw volcano plot
    # Calculate clr
    clr_table <- apply(feature_table, 2, compositions::clr)

    # enrich group
    enrich_group <- function(feature_abd, group) {
      abd_split <- split(feature_abd, group)
      abd_mean_group <- vapply(abd_split, mean, FUN.VALUE = 0.0)
      enrich_group <- names(abd_split)[which.max(abd_mean_group)]

      enrich_group
    }
    group_enriched <- vapply(
      data.frame(t(clr_table)),
      enrich_group,
      FUN.VALUE = character(1),
      group = meta_data %>% dplyr::pull(group_var)
    )

    # Calculate clr mean difference
    eff_size <- apply(clr_table, 1, function(y){
      stats::lm(y ~ x, data = data.frame(y = y,
                                  x = meta_data %>% dplyr::pull(group_var),
                                  check.names = FALSE))$coef[-1]
      }
    )

    if(is.matrix(eff_size)){
      # Data frame for the figure
      dat_fig <- data.frame(FeatureID = out$FeatureID,
                            t(eff_size),
                            y = out$W,
                            check.names = FALSE) %>%
        dplyr::mutate(zero_ind = factor(ifelse(is.infinite(y), "Yes", "No"),
                                        levels = c("Yes", "No"))) %>%
        tidyr::gather(key = group, value = x, rownames(eff_size))
      # Replcace "x" to the name of covariate
      dat_fig$group <- sapply(dat_fig$group, function(x) gsub("x", paste0(group_var, " = "), x))
      # Replace Inf by (n_taxa - 1) for structural zeros
      dat_fig$y <- replace(dat_fig$y, is.infinite(dat_fig$y), n_taxa - 1)

      fig <- ggplot(data = dat_fig, aes(x = x, y = y))+
        geom_point(aes(color = zero_ind))+
        facet_wrap(~ group)+
        labs(x = "CLR mean difference", y = "W statistic")+
        scale_color_discrete(name = "Structural zero", drop = FALSE)+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5),
              legend.position = "top",
              strip.background = element_rect(fill = "white"))
    }else{
      # Data frame for the figure
      dat_fig <- data.frame(taxa_id = out$taxa_id,
                            x = eff_size,
                            y = out$W) %>%
        dplyr::mutate(zero_ind = factor(ifelse(is.infinite(y), "Yes", "No"),
                                        levels = c("Yes", "No")))
      # Replace Inf by (n_taxa - 1) for structural zeros
      dat_fig$y <- replace(dat_fig$y, is.infinite(dat_fig$y), n_taxa - 1)

      fig <- ggplot(data = dat_fig, aes(x = x, y = y))+
        geom_point(aes(color = zero_ind))+
        labs(x = "CLR mean difference", y = "W statistic")+
        scale_color_discrete(name = "Structural zero", drop = FALSE)+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5),
              legend.position = "top")
    }

    res_out <-  eff_size %>% t() %>% data.frame() %>%
      tibble::rownames_to_column("FeatureID") %>%
      dplyr::inner_join(group_enriched %>% data.frame() %>%
                          stats::setNames("Enrichment") %>%
                          tibble::rownames_to_column("FeatureID"),
                        by = "FeatureID") %>%
      dplyr::inner_join(out, by = "FeatureID")

    res <- list(out=res_out, fig=fig)
    return(res)
  }

  metadata <- pData(dataset)
  colnames(metadata)[which(colnames(metadata) == GroupVar)] <- "Group"
  profile <- exprs(dataset)

  # factor group
  metadata$Group <- factor(metadata$Group)
  intersect_sid <- dplyr::intersect(rownames(metadata), colnames(profile))
  # Prepare for input data
  colData <- metadata %>% tibble::rownames_to_column("SampleID") %>%
    dplyr::select(dplyr::all_of(c("SampleID", "Group", AdjVar, RandVar))) %>%
    dplyr::filter(SampleID%in%intersect_sid)
  proData <- profile %>% data.frame() %>%
    dplyr::select(dplyr::all_of(colData$SampleID)) %>%
    as.matrix()

  if(!all(colData$SampleID == colnames(proData))){
    stop("Order of sampleID between colData and proData is wrong please check your data")
  }

  prepro_res <- preprocess_ANCOM(
    proData,
    colData,
    sample_var = "SampleID",
    group_var = GroupVar,
    out_cut = 0.05,
    zero_cut = 0.90,
    lib_cut = 100,
    neg_lb = FALSE)
  res <- ANCOM(prepro_res$feature_table,
               prepro_res$meta_data,
               struc_zero = prepro_res$structure_zeros,
               group_var = GroupVar,
               p_adj_method = "BH",
               alpha = Pvalue,
               adj_formula = AdjVar,
               rand_formula = NULL,
               w_cutoff = Wvalue)

  return(res)
}
