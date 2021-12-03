#' @title Differential Expression Analysis by limma Package
#'
#' @description
#' limma requires the count data (a matrix of integer values) to input. The *Voom* would transform the count data for further data analysis.
#'
#' @details 12/2/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param Expression, ExpressionSet; (Required) ExpressionSet object.
#' @param Group_info, Character; design factor(default: "Group").
#' @param Group_name, Character; (Required) the group for comparison.
#' @param Pvalue, Numeric; significant level(default: 0.05).
#' @param Log2FC, Numeric; log2FoldChange(default: 1).
#'
#' @return
#' a list object:
#'   limma results
#'   significant difference with enriched directors
#'
#' @importFrom dplyr %>% select filter intersect all_of
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom stats setNames model.matrix sd
#' @importFrom edgeR DGEList calcNormFactors
#' @import convert
#' @importFrom limma voom lmFit makeContrasts contrasts.fit eBayes topTable
#'
#' @usage DA_Limma(dataset=ExpressionSet, Group_info="Group", Group_name=c("HC", "AA"), Pvalue=0.05, Log2FC=1)
#' @examples
#'
#' data(ExprSet_species_count)
#'
#' Limma_res <- DA_Limma(dataset=ExprSet_species_count, Group_info="Group", Group_name=c("HC", "AA"), Pvalue=0.05, Log2FC=1)
#' Limma_res$res
#'
DA_Limma <- function(dataset=ExprSet_species_count,
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
  colData <- phen %>% tibble::rownames_to_column("SampleID") %>%
    dplyr::select(dplyr::all_of(c("SampleID", "Group"))) %>%
    dplyr::filter(SampleID%in%intersect_sid) %>%
    dplyr::mutate(Group=factor(Group, levels = Group_name)) %>%
    tibble::column_to_rownames("SampleID")
  countData <- profile %>% data.frame() %>%
    dplyr::select(dplyr::all_of(rownames(colData))) %>%
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

  # DGEList object: counts-> colnames:Samples; rownames:Features
  dge_obj <- edgeR::DGEList(counts = countData)
  dge_obj_factors <- edgeR::calcNormFactors(dge_obj)

  # voom
  dat_voom <- limma::voom(counts = dge_obj_factors,
                          design = design,
                          normalize.method = "quantile",
                          plot = TRUE)

  # Linear model to calculate weighted least squares per gene in each group
  fit <- limma::lmFit(object = dat_voom, design = design)
  # Comparison between two groups
  contr <- paste(rev(Group_name), collapse = "-")
  contr.matrix <- limma::makeContrasts(contrasts = contr,
                                levels = design)
  # Estimate contrast for each gene
  fit.contr <- limma::contrasts.fit(fit, contr.matrix)
  # Empirical Bayes smoothing of standard errors
  fit.ebay <- limma::eBayes(fit.contr)
  # DE result table
  limma_res <- limma::topTable(fit.ebay, coef = contr, n=Inf, sort.by = "P")
  print(head(limma_res))

  # enrichment
  limma_enrich <- limma_res %>% data.frame() %>% arrange(logFC, adj.P.Val)
  if(is.null(Log2FC)){
    Log2FC <- with(limma_enrich,
                   mean(abs(logFC)) + 1.5*stats::sd(abs(logFC)))
    message(paste("Threshold of log2Foldchange [Mean+1.5(SD)] is", Log2FC))
  }else{
    Log2FC <- Log2FC
    message(paste("Threshold of log2Foldchange is", Log2FC))
  }

  limma_enrich[which(limma_enrich$logFC >= Log2FC &
                       limma_enrich$adj.P.Val < Pvalue), "Enrichment"] <- Group_name[2]
  limma_enrich[which(limma_enrich$logFC <= -Log2FC &
                       limma_enrich$adj.P.Val < Pvalue), "Enrichment"] <- Group_name[1]
  limma_enrich[which(abs(limma_enrich$logFC) < Log2FC |
                       limma_enrich$adj.P.Val >= Pvalue), "Enrichment"] <- "Nonsignf"
  print(table(limma_enrich$Enrichment))

  res <- list(res=limma_res,
              enrich=limma_enrich)

  return(res)
}
