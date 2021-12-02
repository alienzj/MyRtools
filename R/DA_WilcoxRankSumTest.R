#' @importFrom dplyr %>% select filter intersect inner_join arrange
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom stats setNames glm
#' @import convert
#'
#' @title Differential Expression Analysis by Wilcox-rank-sum test
#'
#' @description
#' Wilcox-rank-sum test is non-parameter test method, and also use for the data with non-normal distribution(Gauss Distribution).
#'
#' @details 12/2/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param Expression, ExpressionSet; (Required) ExpressionSet object.
#' @param Group_info, Character; design factor(default: "Group").
#' @param Group_name, Character; (Required) the group for comparison.
#' @param Pvalue, Numeric; significant level(default: 0.05).
#'
#' @return
#'   significant difference with enriched directors
#'
#' @usage DA_Wilcox(dataset=ExpressionSet, Group_info="Group", Group_name=c("HC", "AA"), Pvalue=0.05, Log2FC=1)
#' @examples
#'
#' data(ExprSet_species)
#'
#' Wilcox_res <- DA_Wilcox(dataset=ExprSet_species, Group_info="Group", Group_name=c("HC", "AA"), Pvalue=0.05, Log2FC=1)
#' Wilcox_res$res
#'
DA_Wilcox <- function(dataset=ExprSet_species,
                      Group_info="Group",
                      Group_name=c("HC", "AA"),
                      Pvalue=0.05){
  library(dplyr)
  library(convert)
  library(tibble)
  load("data/ExprSet_species.rda")
  dataset=ExprSet_species
  Group_info="Group"
  Group_name=c("HC", "AA")
  Pvalue=0.05

  metadata <- pData(dataset)
  colnames(metadata)[which(colnames(metadata) == Group_info)] <- "Group"
  profile <- exprs(dataset)
  if(!any(profile %% 1 == 0)){
    stop("The input matrix is not integer matrix please Check it")
  }

  # choose group
  phen <- metadata %>% dplyr::filter(Group%in%Group_name)
  intersect_sid <- dplyr::intersect(rownames(phen), colnames(profile))
  # Prepare for input data
  colData <- phen %>% tibble::rownames_to_column("SampleID") %>%
    dplyr::select(SampleID, Group) %>%
    dplyr::filter(SampleID%in%intersect_sid) %>%
    dplyr::mutate(Group=factor(Group, levels = Group_name)) %>%
    tibble::column_to_rownames("SampleID")
  countData <- profile %>% data.frame() %>%
    dplyr::select(all_of(rownames(colData))) %>%
    as.matrix()

  if(!all(rownames(colData) == colnames(countData))){
    stop("Order of sampleID between colData and CountData is wrong please check your data")
  }

  # run wilcox rank sum test
  wilcox_res <- apply(countData, 1, function(x, y){
    dat <- data.frame(value=as.numeric(x), group=y)
    mn <- tapply(dat$value, dat$group, mean) %>%
      data.frame() %>% setNames("value") %>%
      tibble::rownames_to_column("Group")
    mn1 <- with(mn, mn[Group%in%Group_name[1], "value"])
    mn2 <- with(mn, mn[Group%in%Group_name[2], "value"])
    # rank
    rk <- rank(dat$value)
    rnk <- signif(tapply(rk, dat$group, mean), 4)

    Log2FC_rank <- log2(rnk[1]/rnk[2])
    rest <- wilcox.test(data = dat, value ~ group)
    #res <- c(Log2FC, rest$statistic, rest$p.value)
    res <- c(mn1, mn2, rnk, Log2FC_rank, rest$statistic, rest$p.value)
    return(res)
  }, colData$Group) %>%
    t() %>% data.frame() %>%
    stats::setNames(c(paste0("mean_", c(1, 2)), paste0("rank_mean_", c(1, 2)),
               "Log2FC_rank", "Statistic", "P.value")) %>%
    tibble::rownames_to_column("FeatureID") %>%
    mutate(adj.P.Val=p.adjust(as.numeric(P.value), method = "BH"))

  # Number & Block
  dat_status <- table(colData$Group)
  dat_status_number <- as.numeric(dat_status)
  dat_status_name <- names(dat_status)
  wilcox_res$Block <- paste(paste(dat_status_number[1], dat_status_name[1], sep = "_"),
                     "vs",
                     paste(dat_status_number[2], dat_status_name[2], sep = "_"))
  # Enrichment
  wilcox_res$Enrichment_mean <- NA
  for(i in 1:nrow(wilcox_res)){
    if (wilcox_res$mean_1[i] > wilcox_res$mean_2[i] & wilcox_res$adj.P.Val[i] < Pvalue) {
      wilcox_res$Enrichment_mean[i] <- Group_name[1]
    } else if (wilcox_res$mean_1[i] < wilcox_res$mean_2[i] & wilcox_res$adj.P.Val[i] < Pvalue) {
      wilcox_res$Enrichment_mean[i] <- Group_name[2]
    } else if (wilcox_res$adj.P.Val[i] > Pvalue | wilcox_res$mean_1[i] == wilcox_res$mean_2[i]){
      wilcox_res$Enrichment_mean[i] <- "Nonsignif"
    }
  }
  wilcox_res$Enrichment_rank <- NA
  for(i in 1:nrow(wilcox_res)){
    if (wilcox_res$rank_mean_1[i] > wilcox_res$rank_mean_2[i] & wilcox_res$adj.P.Val[i] < Pvalue) {
      wilcox_res$Enrichment_rank[i] <- Group_name[1]
    } else if (wilcox_res$rank_mean_1[i] < wilcox_res$rank_mean_2[i] & wilcox_res$adj.P.Val[i] < Pvalue) {
      wilcox_res$Enrichment_rank[i] <- Group_name[2]
    } else if (wilcox_res$adj.P.Val[i] > Pvalue | wilcox_res$rank_mean_1[i] == wilcox_res$rank_mean_2[i]){
      wilcox_res$Enrichment_rank[i] <- "Nonsignif"
    }
  }

  colnames(wilcox_res)[2:3] <- paste0(Group_name, "_mean")
  colnames(wilcox_res)[4:5] <- paste0(Group_name, "_rank_mean")
  res_temp <- wilcox_res %>% dplyr::select(FeatureID, Block, adj.P.Val,
                                    Enrichment_mean, Enrichment_rank,
                                    everything()) %>% dplyr::arrange(adj.P.Val)

  # glm result for odd ratios 95%CI
  mdat <- dplyr::inner_join(colData %>% tibble::rownames_to_column("SampleID") %>%
                       dplyr::select(all_of(c("SampleID", "Group"))),
                     countData %>% t() %>% data.frame() %>%
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

  res <- inner_join(res_temp, Odd, by="FeatureID")
  return(res)
}