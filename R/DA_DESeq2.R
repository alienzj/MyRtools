#' @title Differential Expression Analysis by DEseq2 Package
#'
#' @description
#' DEseq2 requires the count data (a matrix of integer values) to input. The normalization method is to use the standard factor vector per feature.
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
#'   DESeq object
#'   DESeq results
#'   significant difference with enriched directors
#'
#' @export
#'
#' @importFrom dplyr %>% select filter intersect all_of mutate
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom stats setNames sd
#' @importFrom Biobase pData exprs
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#'
#' @usage DA_DESeq2(dataset=ExpressionSet, Group_info="Group", Group_name=c("HC", "AA"), Pvalue=0.05, Log2FC=1)
#' @examples
#'
#' \donttest{
#' data(ExprSet_species_count)
#'
#' DESeq2_res <- DA_DESeq2(dataset=ExprSet_species_count, Group_info="Group", Group_name=c("HC", "AA"), Pvalue=0.05, Log2FC=1)
#' DESeq2_res$res
#' }
#'
DA_DESeq2 <- function(dataset=ExprSet_species_count,
                      Group_info="Group",
                      Group_name=c("HC", "AA"),
                      Pvalue=0.05,
                      Log2FC=1){

  metadata <- Biobase::pData(dataset)
  colnames(metadata)[which(colnames(metadata) == Group_info)] <- "Group"
  profile <- Biobase::exprs(dataset)
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

  if(!any(rownames(colData) == colnames(countData))){
    stop("Order of sampleID between colData and CountData is wrong please check your data")
  }

  # build Matrix
  ddsm <- DESeq2::DESeqDataSetFromMatrix(
                                 countData=countData,
                                 colData=colData,
                                 design=~ Group)
  # run Deseq2
  dds <- DESeq2::DESeq(ddsm)

  # extract the results
  DESeq_res <- DESeq2::results(dds, contrast = c("Group", rev(Group_name)))
  print(head(DESeq_res))

  # enrichment
  DESeq_enrich <- DESeq_res %>% data.frame() %>% arrange(log2FoldChange, padj)
  if(is.null(Log2FC)){
    Log2FC <- with(DESeq_enrich,
                   mean(abs(log2FoldChange)) + 1.5*stats::sd(abs(log2FoldChange)))
    message(paste("Threshold of log2Foldchange [Mean+1.5(SD)] is", Log2FC))
  }else{
    Log2FC <- Log2FC
    message(paste("Threshold of log2Foldchange is", Log2FC))
  }

  DESeq_enrich[which(DESeq_enrich$log2FoldChange >= Log2FC &
                     DESeq_enrich$padj < Pvalue), "Enrichment"] <- Group_name[2]
  DESeq_enrich[which(DESeq_enrich$log2FoldChange <= -Log2FC &
                     DESeq_enrich$padj < Pvalue), "Enrichment"] <- Group_name[1]
  DESeq_enrich[which(abs(DESeq_enrich$log2FoldChange) < Log2FC |
                     DESeq_enrich$padj >= Pvalue), "Enrichment"] <- "Nonsignf"
  print(table(DESeq_enrich$Enrichment))

  res <- list(dds=dds,
              res=DESeq_res,
              enrich=DESeq_enrich)

  return(res)
}
