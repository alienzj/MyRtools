#' @importFrom dplyr %>% select filter intersect
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom stats setNames model.matrix sd
#' @importFrom limma plotMDS
#' @import convert
#' @importFrom edgeR DGEList calcNormFactors estimateGLMCommonDisp estimateGLMTrendedDisp estimateGLMTagwiseDisp glmFit glmLRT topTags
#'
#' @title Differential Expression Analysis by edgeR Package
#'
#' @description
#' edgeR requires the count data (a matrix of integer values) to input.
#'
#' @details 12/2/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param Expression, ExpressionSet; (Required) ExpressionSet object.
#' @param Group_info, Character; design factor(default: "Group").
#' @param Group_name, Character; the group for comparison.
#' @param Pvalue, Numeric; significant level(default: 0.05).
#' @param Log2FC, Numeric; log2FoldChange(default: 1).
#'
#' @return
#' a list object:
#'   edgeR object
#'   edgeR results
#'   significant difference with enriched directors
#'
#' @usage DA_edgeR(dataset=ExpressionSet, Group_info="Group", Pvalue=0.05, Log2FC=1)
#' @examples
#'
#' data(ExprSet_species_count)
#'
#' Limma_res <- DA_edgeR(dataset=ExprSet_species_count, Group_info="Group", Pvalue=0.05, Log2FC=1)
#' Limma_res$res
#'
DA_edgeR <- function(dataset=ExprSet_species_count,
                     Group_info="Group",
                     Group_name=c("HC", "AA"),
                     Pvalue=0.05,
                     Log2FC=1){

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
  colData <- phen %>% rownames_to_column("SampleID") %>%
    dplyr::select(SampleID, Group) %>%
    dplyr::filter(SampleID%in%intersect_sid) %>%
    dplyr::mutate(Group=factor(Group, levels = Group_name)) %>%
    column_to_rownames("SampleID")
  countData <- profile %>% data.frame() %>%
    dplyr::select(all_of(rownames(colData))) %>%
    as.matrix()
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

  # enrichment
  edgeR_enrich <- edgeR_res %>% data.frame() %>% arrange(logFC, FDR)
  if(is.null(Log2FC)){
    Log2FC <- with(edgeR_enrich,
                   mean(abs(logFC)) + 1.5*stats::sd(abs(logFC)))
    message(paste("Threshold of log2Foldchange [Mean+1.5(SD)] is", Log2FC))
  }else{
    Log2FC <- Log2FC
    message(paste("Threshold of log2Foldchange is", Log2FC))
  }

  edgeR_enrich[which(edgeR_enrich$logFC >= Log2FC &
                       edgeR_enrich$FDR < Pvalue), "Enrichment"] <- Group_name[2]
  edgeR_enrich[which(edgeR_enrich$logFC <= -Log2FC &
                       edgeR_enrich$FDR < Pvalue), "Enrichment"] <- Group_name[1]
  edgeR_enrich[which(abs(edgeR_enrich$logFC) < Log2FC |
                       edgeR_enrich$FDR >= Pvalue), "Enrichment"] <- "Nonsignf"
  print(table(edgeR_enrich$Enrichment))

  res <- list(res=edgeR_res,
              enrich=edgeR_enrich)

  return(res)
}
