#' @title Differential Expression Analysis by T test
#'
#' @description
#' T test is parameter test method, and also use for the data with normal distribution(Gauss Distribution).
#'
#' @details 12/31/2021 Guangzhou China
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
#'   significant difference with enriched directors:
#' * Features
#' * Block
#' * Enrichment
#' * AdjustPvalue
#' * Pvalue
#' * Log2FC
#' * Statistic
#' * Median (All/Groups)
#' * Odds Ratio (95% CI)
#'
#' @export
#'
#' @importFrom dplyr %>% select filter intersect inner_join arrange everything all_of mutate
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom stats setNames glm
#' @importFrom Biobase pData exprs
#'
#' @usage run_Ttest(dataset=ExpressionSet,
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
#' data(ExprSetRawCount)
#'
#' t_res <- run_Ttest(dataset=ExprSetRawCount, Group_info="Group", Group_name=c("HC", "AA"), Pvalue=0.05, Log2FC=1)
#' t_res$res
#' }
#'
run_Ttest <- function(dataset=ExprSetRawCount,
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

  # run test
  t_res <- apply(proData, 1, function(x, y){
    dat <- data.frame(value=as.numeric(x), group=y)
    # median value
    mn <- tapply(dat$value, dat$group, median) %>%
      data.frame() %>% setNames("value") %>%
      tibble::rownames_to_column("Group")
    mn1 <- with(mn, mn[Group%in%Group_name[1], "value"])
    mn2 <- with(mn, mn[Group%in%Group_name[2], "value"])
    mnall <- median(dat$value)
    log2fc <- log2(mn1/mn2)

    rest <- t.test(data = dat, value ~ group, paired = FALSE)
    res <- c(mnall, mn1, mn2, log2fc, rest$statistic, rest$p.value)
    return(res)
  }, colData$Group) %>%
    t() %>% data.frame() %>%
    stats::setNames(c("Median_Abundance", paste0("Median_", c(1, 2)),
                      "Log2FoldChange", "Statistic", "Pvalue")) %>%
    tibble::rownames_to_column("FeatureID") %>%
    mutate(AdjPval=p.adjust(as.numeric(Pvalue), method = "BH"))

  # Number of Group
  dat_status <- table(colData$Group)
  dat_status_number <- as.numeric(dat_status)
  dat_status_name <- names(dat_status)
  t_res$Block <- paste(paste(dat_status_number[1], dat_status_name[1], sep = "_"),
                     "vs",
                     paste(dat_status_number[2], dat_status_name[2], sep = "_"))

  # log2Foldchange exists or not
  index <- which(is.infinite(t_res$Log2FoldChange) | is.nan(t_res$Log2FoldChange))
  t_res_Inf <- t_res[index, ]
  t_res_NoInf <- t_res[-index, ]

  # Enrichment
  if(nrow(t_res_Inf) != 0){
    t_res_Inf[which(t_res_Inf$Median_1 > t_res_Inf$Median_2 & t_res_Inf$AdjPval < Pvalue), "Enrichment"] <- Group_name[1]
    t_res_Inf[which(t_res_Inf$Median_1 < t_res_Inf$Median_2 & t_res_Inf$AdjPval < Pvalue), "Enrichment"] <- Group_name[2]
    t_res_Inf[which(t_res_Inf$Median_1 == t_res_Inf$Median_2 | t_res_Inf$AdjPval >= Pvalue), "Enrichment"] <- "Nonsignif"
  }

  if(nrow(t_res_NoInf) != 0){
    t_res_NoInf[which(t_res_NoInf$Log2FoldChange > Log2FC & t_res_NoInf$AdjPval < Pvalue), "Enrichment"] <- Group_name[1]
    t_res_NoInf[which(t_res_NoInf$Log2FoldChange < -Log2FC & t_res_NoInf$AdjPval < Pvalue), "Enrichment"] <- Group_name[2]
    t_res_NoInf[which(abs(t_res_NoInf$Log2FoldChange) <= Log2FC | t_res_NoInf$AdjPval >= Pvalue), "Enrichment"] <- "Nonsignif"
  }

  t_res_merge <- rbind(t_res_Inf, t_res_NoInf)

  # Rename
  colnames(t_res_merge)[2:4] <- c("Median Abundance\n(All)", paste0("Median Abundance\n", Group_name))
  t_res_temp <- t_res_merge %>% dplyr::select(FeatureID, Block, Enrichment,
                                           AdjPval, Pvalue, Log2FoldChange, Statistic,
                                    dplyr::everything()) %>% dplyr::arrange(AdjPval)

  # 95% CI Odd Ratio
  res_odd <- run_OddRatio(colData, proData, Group_name)

  # Merge
  res <- dplyr::inner_join(t_res_temp, res_odd, by="FeatureID")

  return(res)
}


#' @title Differential Expression Analysis by Paired T test
#'
#' @description
#' Paired T test is parameter test method, and also use for the data with normal distribution(Gauss Distribution).
#'
#' @details 12/31/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param Expression, ExpressionSet; (Required) ExpressionSet object.
#' @param trim, Character; filter to apply.(default: trim="none").
#' @param transform, Character; transformation to apply.(default: tranform="none").
#' @param normalize, Character; normalization to apply.(default: normalize="none").
#' @param Group_info, Character; design factor(default: "Group").
#' @param Group_name, Character; (Required) the group for comparison.
#' @param Pair_ID, Character; (Required) the paired ID.
#' @param Pvalue, Numeric; significant level(default: 0.05).
#' @param Log2FC, Numeric; log2FoldChange(default: 1).
#'
#' @return
#'   significant difference with enriched directors:
#' * Features
#' * Block
#' * Enrichment
#' * AdjustPvalue
#' * Pvalue
#' * Log2FC
#' * Statistic
#' * Median (All/Groups)
#' * Odds Ratio (95% CI)
#'
#' @export
#'
#' @importFrom dplyr %>% select filter intersect inner_join arrange everything all_of mutate
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom stats setNames glm
#' @importFrom Biobase pData exprs
#'
#' @usage run_Ttest_Paired(dataset=ExpressionSet,
#'                   trim="none",
#'                   transform="none",
#'                   normalize="none",
#'                   Group_info="Group",
#'                   Group_name=c("HC", "AA"),
#'                   Pair_ID="PID",
#'                   Pvalue=0.05,
#'                   Log2FC=1)
#' @examples
#'
#' \donttest{
#' data(ExprSetRawCount)
#'
#' t_res <- run_Ttest_Paired(dataset=ExprSetRawCount, Group_info="Group", Group_name=c("HC", "AA"), Pair_ID="PID", Pvalue=0.05, Log2FC=1)
#' t_res$res
#' }
#'
run_Ttest_Paired <- function(dataset=ExprSetRawCount,
                              trim="none",
                              transform="none",
                              normalize="none",
                              Group_info="Group",
                              Group_name=c("HC", "AA"),
                              Pair_ID="PID",
                              Pvalue=0.05,
                              Log2FC=1){

  # preprocess
  dataset_processed <- get_processedExprSet(dataset=dataset,
                                            trim=trim,
                                            transform=transform,
                                            normalize=normalize)

  metadata <- Biobase::pData(dataset_processed)
  colnames(metadata)[which(colnames(metadata) == Group_info)] <- "Group"
  colnames(metadata)[which(colnames(metadata) == Pair_ID)] <- "PID"
  profile <- Biobase::exprs(dataset_processed)

  # Choose group
  phen <- metadata %>% dplyr::filter(Group%in%Group_name)
  if(length(unique(phen$Group)) != 2){
    stop("Levels of Group must be 2")
  }
  # choose PID
  intersect_pid <- dplyr::intersect(unique(phen$PID[which(phen$Group == Group_name[1])]),
                                    unique(phen$PID[which(phen$Group == Group_name[2])]))
  if(length(intersect_pid) == 0){
    stop("There are no common PID Please check your metadata")
  }else{
    number_pid <- length(intersect_pid)
    message(paste0("The number of Paired Subjects is ", number_pid))
  }
  phen_paired <- phen %>% dplyr::filter(PID%in%intersect_pid)

  intersect_sid <- dplyr::intersect(rownames(phen_paired), colnames(profile))
  # Prepare for input data
  colData <- phen_paired %>% tibble::rownames_to_column("SampleID") %>%
    dplyr::select(dplyr::all_of(c("PID", "SampleID", "Group"))) %>%
    dplyr::filter(SampleID%in%intersect_sid) %>%
    dplyr::mutate(Group=factor(Group, levels = Group_name)) %>%
    dplyr::arrange(PID, SampleID) %>%
    tibble::column_to_rownames("SampleID")
  proData <- profile %>% data.frame() %>%
    dplyr::select(dplyr::all_of(rownames(colData))) %>%
    as.matrix()

  if(!all(rownames(colData) == colnames(proData))){
    stop("Order of sampleID between colData and proData is wrong please check your data")
  }

  # run wilcox rank sum test
  t_res <- apply(proData, 1, function(x, y){
    dat <- data.frame(value=as.numeric(x), group=y$Group, pid=y$PID) %>%
      dplyr::arrange(pid, group)
    # median value
    mn <- tapply(dat$value, dat$group, median) %>%
      data.frame() %>% setNames("value") %>%
      tibble::rownames_to_column("Group")
    mn1 <- with(mn, mn[Group%in%Group_name[1], "value"])
    mn2 <- with(mn, mn[Group%in%Group_name[2], "value"])
    mnall <- median(dat$value)
    log2fc <- log2(mn1/mn2)

    rest <- t.test(data = dat, value ~ group, paired = TRUE)
    res <- c(mnall, mn1, mn2, log2fc, rest$statistic, rest$p.value)
    return(res)
  }, colData) %>%
    t() %>% data.frame() %>%
    stats::setNames(c("Median_Abundance", paste0("Median_", c(1, 2)),
                      "Log2FoldChange", "Statistic", "Pvalue")) %>%
    tibble::rownames_to_column("FeatureID") %>%
    mutate(AdjPval=p.adjust(as.numeric(Pvalue), method = "BH"))

  # Number of Group
  dat_status <- table(colData$Group)
  dat_status_number <- as.numeric(dat_status)
  dat_status_name <- names(dat_status)
  t_res$Block <- paste("Paired", paste(dat_status_number[1], dat_status_name[1], sep = "_"),
                       "vs",
                       paste(dat_status_number[2], dat_status_name[2], sep = "_"))

  # log2Foldchange exists or not
  index <- which(is.infinite(t_res$Log2FoldChange) | is.nan(t_res$Log2FoldChange))
  t_res_Inf <- t_res[index, ]
  t_res_NoInf <- t_res[-index, ]

  # Enrichment
  if(nrow(t_res_Inf) != 0){
    t_res_Inf[which(t_res_Inf$Median_1 > t_res_Inf$Median_2 & t_res_Inf$AdjPval < Pvalue), "Enrichment"] <- Group_name[1]
    t_res_Inf[which(t_res_Inf$Median_1 < t_res_Inf$Median_2 & t_res_Inf$AdjPval < Pvalue), "Enrichment"] <- Group_name[2]
    t_res_Inf[which(t_res_Inf$Median_1 == t_res_Inf$Median_2 | t_res_Inf$AdjPval >= Pvalue), "Enrichment"] <- "Nonsignif"
  }

  if(nrow(t_res_NoInf) != 0){
    t_res_NoInf[which(t_res_NoInf$Log2FoldChange > Log2FC & t_res_NoInf$AdjPval < Pvalue), "Enrichment"] <- Group_name[1]
    t_res_NoInf[which(t_res_NoInf$Log2FoldChange < -Log2FC & t_res_NoInf$AdjPval < Pvalue), "Enrichment"] <- Group_name[2]
    t_res_NoInf[which(abs(t_res_NoInf$Log2FoldChange) <= Log2FC | t_res_NoInf$AdjPval >= Pvalue), "Enrichment"] <- "Nonsignif"
  }

  t_res_merge <- rbind(t_res_Inf, t_res_NoInf)

  # Rename
  colnames(t_res_merge)[2:4] <- c("Median Abundance\n(All)", paste0("Median Abundance\n", Group_name))
  t_res_temp <- t_res_merge %>% dplyr::select(FeatureID, Block, Enrichment,
                                        AdjPval, Pvalue, Log2FoldChange, Statistic,
                                        dplyr::everything()) %>% dplyr::arrange(AdjPval)

  # 95% CI Odd Ratio
  res_odd <- run_OddRatio(colData, proData, Group_name)

  # Merge
  res <- dplyr::inner_join(t_res_temp, res_odd, by="FeatureID")

  return(res)
}

