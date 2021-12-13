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
#' @param trim, Character; filter to apply.(default: trim="none").
#' @param transform, Character; transformation to apply.(default: tranform="none").
#' @param normalize, Character; normalization to apply.(default: normalize="none").
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
#' @usage run_Wilcox(dataset=ExpressionSet,
#'                   trim="none",
#'                   transform="none",
#'                   normalize="none",
#'                   Group_info="Group",
#'                   Group_name=c("HC", "AA"),
#'                   Pvalue=0.05,
#'                   Log2FC=1)
#' @examples
#'
#' \donttest{
#' data(ExprSetRawRB)
#'
#' TWilcox_res <- run_TWilcox(dataset=ExprSetRawRB, Group_info="Group", Group_name=c("HC", "AA"), Pvalue=0.05, Log2FC=1)
#' TWilcox_res$res
#' }
#'
run_TWilcox <- function(dataset=ExprSetRawCount,
                        trim="none",
                        transform="none",
                        normalize="none",
                        Group_info="Group",
                        Group_name=c("HC", "AA"),
                        Pvalue=0.05,
                        Log2FC=1){

  # preprocess
  dataset_processed <- get_processedExprSet(dataset=dataset,
                                            trim=trim,
                                            transform=transform,
                                            normalize=normalize)
  metadata <- Biobase::pData(dataset_processed)
  colnames(metadata)[which(colnames(metadata) == Group_info)] <- "Group"
  profile <- Biobase::exprs(dataset_processed)

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

      # median(Original Value)
      mn <- tapply(dat$value, dat$group, median) %>%
        data.frame() %>% setNames("value") %>%
        tibble::rownames_to_column("Group")
      mn1 <- with(mn, mn[Group%in%Group_name[1], "value"])
      mn2 <- with(mn, mn[Group%in%Group_name[2], "value"])
      mnall <- median(dat$value)

      # Fold2Change: geometricmean
      dat$value_scale <- scale(dat$value, center = TRUE, scale = TRUE)
      mn_GM <- tapply(dat$value_scale, dat$group, compositions::geometricmean) %>%
        data.frame() %>% stats::setNames("value") %>%
        rownames_to_column("Group")
      mn_GM1 <- with(mn_GM, mn_GM[Group%in%Group_name[1], "value"])
      mn_GM2 <- with(mn_GM, mn_GM[Group%in%Group_name[2], "value"])
      Log2FC_GM <- log2(mn_GM1/mn_GM2)

      # pvalue
      rest <- rstatix::t_test(data = dat, value ~ group)

      return(c(mnall, mn1, mn2,
               mn_GM1, mn_GM2,Log2FC_GM,
               rest$statistic, rest$p))
    }, colData$Group) %>%
      t() %>% data.frame() %>%
    stats::setNames(c("Median_Abundance", paste0("Median_", c(1, 2)),
                      paste0("geometricmean_", c(1, 2)),
                      "Log2FoldChange", "Statistic", "Pvalue"))

    Normal_res <- Welch_res %>%
      # dplyr::filter(!is.nan(Log2FoldChange)) %>%
      # dplyr::filter(!is.infinite(Log2FoldChange)) %>%
      tibble::rownames_to_column("FeatureID") %>%
      dplyr::arrange(desc(abs(Log2FoldChange)), Pvalue)
  }else{
    Normal_res <- data.frame()
  }

  # Non-Normally Distributed Features
  if(nrow(Non_Normal_prof) != 0){
    Wilcox_res <- apply(Non_Normal_prof, 1, function(x, y){
      dat <- data.frame(value=as.numeric(x), group=y)
      # median(Original Value)
      mn <- tapply(dat$value, dat$group, median) %>%
        data.frame() %>% setNames("value") %>%
        tibble::rownames_to_column("Group")
      mn1 <- with(mn, mn[Group%in%Group_name[1], "value"])
      mn2 <- with(mn, mn[Group%in%Group_name[2], "value"])
      mnall <- median(dat$value)

      # Fold2Change: geometricmean
      dat$value_scale <- scale(dat$value, center = TRUE, scale = TRUE)
      mn_GM <- tapply(dat$value_scale, dat$group, compositions::geometricmean) %>%
        data.frame() %>% stats::setNames("value") %>%
        rownames_to_column("Group")
      mn_GM1 <- with(mn_GM, mn_GM[Group%in%Group_name[1], "value"])
      mn_GM2 <- with(mn_GM, mn_GM[Group%in%Group_name[2], "value"])
      Log2FC_GM <- log2(mn_GM1/mn_GM2)

      # pvalue
      rest <- stats::wilcox.test(data = dat, value ~ group)

      return(c(mnall, mn1, mn2,
               mn_GM1, mn_GM2,Log2FC_GM,
               rest$statistic, rest$p.value))
    }, colData$Group) %>%
      t() %>% data.frame() %>%
      stats::setNames(c("Median_Abundance", paste0("Median_", c(1, 2)),
                        paste0("geometricmean_", c(1, 2)),
                        "Log2FoldChange", "Statistic", "Pvalue"))

    Non_Normal_res <- Wilcox_res %>%
      # dplyr::filter(!is.nan(Log2FoldChange)) %>%
      # dplyr::filter(!is.infinite(Log2FoldChange)) %>%
      tibble::rownames_to_column("FeatureID") %>%
      dplyr::arrange(desc(abs(Log2FoldChange)), Pvalue)
  }else{
    Non_Normal_res <- data.frame()
  }

  # Number & Block
  t_res <- rbind(Normal_res, Non_Normal_res)
  t_res$AdjPval <- p.adjust(as.numeric(t_res$Pvalue), method = "BH")

  dat_status <- table(colData$Group)
  dat_status_number <- as.numeric(dat_status)
  dat_status_name <- names(dat_status)
  t_res$Block <- paste(paste(dat_status_number[1], dat_status_name[1], sep = "_"),
                     "vs",
                     paste(dat_status_number[2], dat_status_name[2], sep = "_"))

  # Enrichment geometricmean
  t_res[which(t_res$Log2FoldChange > Log2FC & t_res$AdjPval < Pvalue), "Enrichment"] <- Group_name[1]
  t_res[which(t_res$Log2FoldChange < -Log2FC & t_res$AdjPval < Pvalue), "Enrichment"] <- Group_name[2]
  t_res[which(abs(t_res$Log2FoldChange) <= Log2FC | t_res$AdjPval >= Pvalue), "Enrichment"] <- "Nonsignif"

  # Rename
  colnames(t_res)[2:4] <- c("Median Abundance\n(All)", paste0("Median Abundance\n", Group_name))
  colnames(t_res)[5:6] <- paste0("Geometricmean Abundance\n", Group_name)
  t_res_temp <- t_res %>% dplyr::select(FeatureID, Block, Enrichment,
                                        AdjPval, Pvalue, Log2FoldChange, Statistic,
                                        dplyr::everything()) %>% dplyr::arrange(AdjPval)

  # 95% CI Odd Ratio
  res_odd <- run_OddRatio(colData, proData, Group_name)

  # Merge
  res <- dplyr::inner_join(t_res_temp, res_odd, by="FeatureID")

  return(res)
}
