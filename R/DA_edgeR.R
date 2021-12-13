#' @title Differential Expression Analysis by edgeR Package
#'
#' @description
#' edgeR requires the count data (a matrix of integer values) to input.
#'
#' @details 12/2/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param Expression, ExpressionSet; (Required) ExpressionSet object.
#' @param trim, Character; filter to apply.(default: trim="none").
#' @param Group_info, Character; design factor(default: "Group").
#' @param Group_name, Character; (Required) the group for comparison.
#' @param Pvalue, Numeric; significant level(default: 0.05).
#' @param Log2FC, Numeric; log2FoldChange(default: 1).
#'
#' @return
#' a list object:
#'   edgeR results
#'   significant difference with enriched directors
#'
#' @export
#'
#' @importFrom dplyr %>% select filter intersect all_of mutate
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom stats setNames model.matrix sd
#' @importFrom limma plotMDS
#' @importFrom Biobase pData exprs
#' @importFrom edgeR DGEList calcNormFactors estimateGLMCommonDisp estimateGLMTrendedDisp estimateGLMTagwiseDisp glmFit glmLRT topTags
#'
#' @usage run_edgeR(dataset=ExpressionSet,
#'                  trim="none",
#'                  Group_info="Group",
#'                  Group_name=c("HC", "AA"),
#'                  Pvalue=0.05,
#'                  Log2FC=1)
#' @examples
#'
#' \donttest{
#' data(ExprSetRawCount)
#'
#' edgeR_res <- run_edgeR(dataset=ExprSetRawCount, Group_info="Group", Group_name=c("HC", "AA"), Pvalue=0.05, Log2FC=1)
#' edgeR_res$res
#' }
#'
run_edgeR <- function(dataset=ExprSetRawCount,
                      trim="none",
                      Group_info="Group",
                      Group_name=c("HC", "AA"),
                      Pvalue=0.05,
                      Log2FC=1){

  # preprocess
  dataset_processed <- get_processedExprSet(dataset=dataset, trim=trim)

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
  countData <- profile %>% data.frame() %>%
    dplyr::select(dplyr::all_of(rownames(colData))) %>%
    as.matrix()

  # Median abundance
  median_res <- apply(countData, 1, function(x, y){
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

  # No zero value for Log transform
  if(any(countData == 0)){
    countData <- countData+1
  }else{
    countData <- countData
  }

  if(!all(rownames(colData) == colnames(countData))){
    stop("Order of sampleID between colData and CountData is wrong please check your data")
  }

  # group Matrix
  design <- stats::model.matrix(~0 + Group, data = colData)
  colnames(design) <- levels(colData$Group)
  rownames(design) <- colnames(countData)

  # Filter data
  keep <- rowSums(edgeR::cpm(countData) > 100) >= 2
  countData <- countData[keep, ]

  # DGEList object: counts-> colnames:Samples; rownames:Features
  dge_obj <- edgeR::DGEList(counts = countData)
  dge_obj_factors <- edgeR::calcNormFactors(dge_obj)
  limma::plotMDS(dge_obj_factors, method="bcv", col=as.numeric(colData$Group))
  legend("bottomleft", Group_name, col=1:2, pch=20)

  # GLM estimates of dispersion
  dge <- edgeR::estimateGLMCommonDisp(dge_obj_factors, design)
  dge <- edgeR::estimateGLMTrendedDisp(dge, design, method="power")
  # You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
  # The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
  dge <- edgeR::estimateGLMTagwiseDisp(dge, design)
  #edgeR::plotBCV(dge)

  # GLM testing for differential expression:
  fit <- edgeR::glmFit(dge, design)
  fit2 <- edgeR::glmLRT(fit, contrast=c(-1, 1)) # Normal:-1; Tumor:1
  edgeR_res <- edgeR::topTags(fit2, n=nrow(countData)) %>% data.frame()
  print(head(edgeR_res))

  colnames(edgeR_res)[c(1:2, 4:5)] <- c("Log2FoldChange", "logCPM", "Pvalue", "AdjPval")
  t_res <- edgeR_res %>% tibble::rownames_to_column("FeatureID") %>%
    dplyr::select(dplyr::all_of(c("FeatureID", "Log2FoldChange",
                                  "logCPM", "Pvalue", "AdjPval"))) %>%
    dplyr::inner_join(median_res, by = "FeatureID")

  # enrichment
  if(is.null(Log2FC)){
    Log2FC <- with(t_res,
                   mean(abs(Log2FoldChange)) + 1.5*stats::sd(abs(Log2FoldChange)))
    message(paste("Threshold of log2Foldchange [Mean+1.5(SD)] is", Log2FC))
  }else{
    Log2FC <- Log2FC
    message(paste("Threshold of log2Foldchange is", Log2FC))
  }
  t_res[which(t_res$Log2FoldChange > Log2FC & t_res$AdjPval < Pvalue), "Enrichment"] <- Group_name[2]
  t_res[which(t_res$Log2FoldChange < -Log2FC & t_res$AdjPval < Pvalue), "Enrichment"] <- Group_name[1]
  t_res[which(abs(t_res$Log2FoldChange) <= Log2FC | t_res$AdjPval >= Pvalue), "Enrichment"] <- "Nonsignif"
  print(table(t_res$Enrichment))

  # Number of Group
  dat_status <- table(colData$Group)
  dat_status_number <- as.numeric(dat_status)
  dat_status_name <- names(dat_status)
  t_res$Block <- paste(paste(dat_status_number[1], dat_status_name[1], sep = "_"),
                       "vs",
                       paste(dat_status_number[2], dat_status_name[2], sep = "_"))
  t_res_temp <- t_res %>% dplyr::select(FeatureID, Block, Enrichment,
                                        AdjPval, Pvalue, Log2FoldChange, logCPM,
                                        dplyr::everything()) %>% dplyr::arrange(AdjPval)
  # 95% CI Odd Ratio
  res_odd <- run_OddRatio(colData, countData, Group_name)

  res <- dplyr::inner_join(t_res_temp, res_odd, by="FeatureID")

  return(res)
}
