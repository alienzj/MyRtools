#' @title Differential Expression Analysis by T test or Wilcox-rank-sum test
#'
#' @description
#' The data of features maybe either non-normal distribution or normal distribution(Gauss Distribution).
#' shapiro.test is used to check the distribution. T-test is for normal distribution, while Wilcox-rank-sum test is for non-normal distribution.
#' The log2FoldChange is calculate by the geometricmean or median of two groups after a linear scale (Zscore).
#'
#' @details 12/3/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param Expression, ExpressionSet; (Required) ExpressionSet object.
#' @param Group_info, Character; design factor(default: "Group").
#' @param Group_name, Character; (Required) the group for comparison.
#' @param Pvalue, Numeric; significant level(default: 0.05).
#' @param Log2FC, Numeric; log2FoldChange(default: 1).
#'
#' @return
#'   significant difference with enriched directors
#'
#' @export
#'
#' @importFrom dplyr %>% select filter intersect inner_join arrange everything all_of mutate
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom stats setNames glm shapiro.test median wilcox.test
#' @importFrom compositions geometricmean
#' @importFrom rstatix t_test
#' @importFrom Biobase pData exprs
#'
#' @usage DA_TWilcox(dataset=ExpressionSet, Group_info="Group", Group_name=c("HC", "AA"), Pvalue=0.05, Log2FC=1)
#' @examples
#'
#' \donttest{
#' data(ExprSet_species)
#'
#' TWilcox_res <- DA_TWilcox(dataset=ExprSet_species, Group_info="Group", Group_name=c("HC", "AA"), Pvalue=0.05, Log2FC=1)
#' TWilcox_res$res
#' }
#'
DA_TWilcox <- function(dataset=ExprSet_species,
                       Group_info="Group",
                       Group_name=c("HC", "AA"),
                       Pvalue=0.05,
                       Log2FC=1){

  metadata <- Biobase::pData(dataset)
  colnames(metadata)[which(colnames(metadata) == Group_info)] <- "Group"
  profile <- Biobase::exprs(dataset)
  # if(!any(profile %% 1 == 0)){
  #   stop("The input matrix is not integer matrix please Check it")
  # }

  # Choose group
  phen <- metadata %>% dplyr::filter(Group%in%Group_name)
  if(length(unique(phen$Group)) != 2){
    stop("Levels of Group must be 2")
  }
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

  # Checking whether the two groups that where tested were normally distributed using the Shapiro-Wilk test
  shapiro_res <- apply(proData, 1, function(x, y){
    dat <- data.frame(value=x, group=y)
    t_res <- tapply(dat$value, dat$group, stats::shapiro.test)
    if(t_res[[1]]$p.value < 0.05){
      res <- FALSE
    }else{
      res <- TRUE
    }
    return(res)
  }, colData$Group) %>%
    data.frame() %>%
    stats::setNames("Normal") %>%
    tibble::rownames_to_column("FeatureID")

  Normal_prof <- proData[shapiro_res$Normal, ]
  Non_Normal_prof <- proData[!shapiro_res$Normal, ]

  # Normally Distributed Features
  if(nrow(Normal_prof) != 0){
    Welch_res <- apply(Normal_prof, 1, function(x, y){
      dat <- data.frame(value=as.numeric(x), group=y)
      # Fold2Change: median
      dat$value_scale <- scale(dat$value, center = TRUE, scale = TRUE)
      mn_median <- tapply(dat$value_scale, dat$group, stats::median) %>%
        data.frame() %>% stats::setNames("value") %>%
        rownames_to_column("Group")
      mn_median1 <- with(mn_median, mn_median[Group%in%group_name[1], "value"])
      mn_median2 <- with(mn_median, mn_median[Group%in%group_name[2], "value"])
      if(all(mn_median1 > 0, mn_median2 > 0)){
        Log2FC_median <- log2(mn_median1/mn_median2)
      }else{
        Log2FC_median <- -log2(mn_median1/mn_median2)
      }
      # Fold2Change: geometricmean
      mn_GM <- tapply(dat$value_scale, dat$group, compositions::geometricmean) %>%
        data.frame() %>% stats::setNames("value") %>%
        rownames_to_column("Group")
      mn_GM1 <- with(mn_GM, mn_GM[Group%in%group_name[1], "value"])
      mn_GM2 <- with(mn_GM, mn_GM[Group%in%group_name[2], "value"])
      Log2FC_GM <- log2(mn_GM1/mn_GM2)

      # pvalue
      rest <- rstatix::t_test(data = dat, value ~ group)

      return(c(Log2FC_median, mn_median1, mn_median2,
               Log2FC_GM, mn_GM1, mn_GM2,
               rest$statistic, rest$p))
    }, colData$Group) %>%
      t() %>% data.frame() %>%
      stats::setNames(c("log2FC_median", paste0("median_", Group_name),
                 "Log2FC_GM", paste0("geometricmean_", Group_name),
                 "Statistic", "P.value"))
    Normal_res <- Welch_res %>%
      # dplyr::filter(!is.nan(Log2FC_GM)) %>%
      # dplyr::filter(!is.infinite(Log2FC_GM)) %>%
      tibble::rownames_to_column("FeatureID") %>%
      dplyr::arrange(desc(abs(Log2FC_GM)), P.value)
  }else{
    Normal_res <- data.frame()
  }

  # Non-Normally Distributed Features
  if(nrow(Non_Normal_prof) != 0){
    Wilcox_res <- apply(Non_Normal_prof, 1, function(x, y){
      dat <- data.frame(value=as.numeric(x), group=y)
      # Fold2Change: median
      dat$value_scale <- scale(dat$value, center = TRUE, scale = TRUE)
      mn_median <- tapply(dat$value_scale, dat$group, stats::median) %>%
        data.frame() %>% stats::setNames("value") %>%
        tibble::rownames_to_column("Group")
      mn_median1 <- with(mn_median, mn_median[Group%in%Group_name[1], "value"])
      mn_median2 <- with(mn_median, mn_median[Group%in%Group_name[2], "value"])
      if(all(mn_median1 > 0, mn_median2 > 0)){
        Log2FC_median <- log2(mn_median1/mn_median2)
      }else{
        Log2FC_median <- -log2(mn_median1/mn_median2)
      }
      # Fold2Change: geometricmean
      mn_GM <- tapply(dat$value_scale, dat$group, compositions::geometricmean) %>%
        data.frame() %>% stats::setNames("value") %>%
        tibble::rownames_to_column("Group")
      mn_GM1 <- with(mn_GM, mn_GM[Group%in%Group_name[1], "value"])
      mn_GM2 <- with(mn_GM, mn_GM[Group%in%Group_name[2], "value"])
      Log2FC_GM <- log2(mn_GM1/mn_GM2)

      # pvalue
      rest <- stats::wilcox.test(data = dat, value ~ group)

      return(c(Log2FC_median, mn_median1, mn_median2,
               Log2FC_GM, mn_GM1, mn_GM2,
               rest$statistic, rest$p.value))
    }, colData$Group) %>%
      t() %>% data.frame() %>%
      stats::setNames(c("log2FC_median", paste0("median_", Group_name),
                 "log2FC_GM", paste0("geometricmean_", Group_name),
                 "Statistic", "P.value"))
    Non_Normal_res <- Wilcox_res %>%
      # dplyr::filter(!is.nan(log2FC_GM)) %>%
      # dplyr::filter(!is.infinite(log2FC_GM)) %>%
      tibble::rownames_to_column("FeatureID") %>%
      dplyr::arrange(desc(abs(log2FC_GM)), P.value)
  }else{
    Non_Normal_res <- data.frame()
  }

  # Number & Block
  all_res <- rbind(Normal_res, Non_Normal_res)
  all_res$adj.P.Val <- p.adjust(as.numeric(all_res$P.value), method = "BH")

  dat_status <- table(colData$Group)
  dat_status_number <- as.numeric(dat_status)
  dat_status_name <- names(dat_status)
  all_res$Block <- paste(paste(dat_status_number[1], dat_status_name[1], sep = "_"),
                     "vs",
                     paste(dat_status_number[2], dat_status_name[2], sep = "_"))
  # Enrichment Meidan
  all_res[which(all_res$log2FC_median > Log2FC & all_res$adj.P.Val < Pvalue), "Enrichment_median"] <- Group_name[1]
  all_res[which(all_res$log2FC_median < -Log2FC & all_res$adj.P.Val < Pvalue), "Enrichment_median"] <- Group_name[2]
  all_res[which(abs(all_res$log2FC_median) <= Log2FC | all_res$adj.P.Val >= Pvalue), "Enrichment_median"] <- "Nonsignif"

  # Enrichment geometricmean
  all_res[which(all_res$log2FC_GM > Log2FC & all_res$adj.P.Val < Pvalue), "Enrichment_GM"] <- Group_name[1]
  all_res[which(all_res$log2FC_GM < -Log2FC & all_res$adj.P.Val < Pvalue), "Enrichment_GM"] <- Group_name[2]
  all_res[which(abs(all_res$log2FC_GM) <= Log2FC | all_res$adj.P.Val >= Pvalue), "Enrichment_GM"] <- "Nonsignif"

  res_temp <- all_res %>% dplyr::select(FeatureID, Block, adj.P.Val, P.value, dplyr::everything()) %>%
    arrange(adj.P.Val)

  # glm result for odd ratios 95%CI
  mdat <- dplyr::inner_join(colData %>% tibble::rownames_to_column("SampleID") %>%
                       dplyr::select(dplyr::all_of(c("SampleID", "Group"))),
                     proData %>% t() %>% data.frame() %>%
                       tibble::rownames_to_column("SampleID"),
                     by = "SampleID") %>%
    tibble::column_to_rownames("SampleID")

  dat_phe <- mdat %>% dplyr::select(Group) %>%
    mutate(Group=ifelse(Group==Group_name[2], 1, 0))
  dat_prf <- mdat %>% dplyr::select(-Group)

  glmFun <- function(GroupN, MarkerN){

    MarkerN[MarkerN==0] <- min(MarkerN[MarkerN!=0])
    dat_glm <- data.frame(group=GroupN,
                          marker=scale(MarkerN, center=TRUE, scale=TRUE)) %>%
      na.omit()
    model <- summary(stats::glm(group ~ marker, data = dat_glm,
                         family = binomial(link = "logit")))
    res <- signif(exp(model$coefficients["marker", 1]) +
                    qnorm(c(0.025,0.5,0.975)) * model$coefficients["marker",1], 2)

    return(res)
  }

  glm_res <- t(apply(dat_prf, 2, function(x, group){
    res <- glmFun(group, as.numeric(x))
    return(res)
  }, dat_phe$Group))

  Odd <- glm_res %>% data.frame() %>%
    setNames(c("upper", "expected","lower")) %>%
    mutate("Odds Ratio (95% CI)" = paste0(expected, " (", lower, ";", upper, ")"))
  Odd$FeatureID <- rownames(glm_res)

  res <- inner_join(res_temp, Odd[, c(4:5)], by="FeatureID")
  return(res)
}
