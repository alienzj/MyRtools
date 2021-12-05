#' @title Association Analysis by Spearman correlation coefficient
#'
#' @description
#'  The correlation among features is calculated by Spearman.
#'
#' @details 12/3/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param Expression, ExpressionSet; (Required) ExpressionSet object.
#' @param Group_info, Character; design factor(default: "Group").
#' @param Group_name, Character; the group for filtering(default: NULL).
#' @param Object, Character; (Required) the variables for correlation(Feature;Meta;Mix).
#'
#' @return
#'  a list of results
#'   Spearman Correlation Coefficient matrix
#'   Spearman Correlation Coefficient table
#'
#' @export
#'
#' @importFrom dplyr %>% select filter intersect inner_join all_of
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom stats na.omit
#' @importFrom Hmisc rcorr
#' @importFrom Biobase pData exprs
#'
#' @usage AA_SimpleCC(dataset=ExpressionSet, Group_info="Group", Group_name=NULL, Object="Feature")
#' @examples
#'
#' \donttest{
#' data(ExprSet_species)
#'
#' SimpleCC_res <- AA_SimpleCC(dataset=ExprSet_species, Group_info="Group", Group_name=NULL, Object="Feature")
#'
#'}
#'
AA_SimpleCC <- function(dataset=ExprSet_species,
                        Group_info="Group",
                        Group_name=NULL,
                        Object="Feature"){

  metadata <- Biobase::pData(dataset)
  colnames(metadata)[which(colnames(metadata) == Group_info)] <- "Group"
  profile <- Biobase::exprs(dataset)

  # Prepare for input data
  if(is.null(Group_name)){
    colData <- metadata
  }else{
    phen <- metadata %>% dplyr::filter(Group%in%Group_name)
    intersect_sid <- dplyr::intersect(rownames(phen), colnames(profile))
    colData <- phen %>% tibble::rownames_to_column("SampleID") %>%
      dplyr::filter(SampleID%in%intersect_sid) %>%
      tibble::column_to_rownames("SampleID")
  }

  proData <- profile %>% data.frame() %>%
    dplyr::select(dplyr::all_of(rownames(colData))) %>%
    as.matrix()

  if(!all(rownames(colData) == colnames(proData))){
    stop("Order of sampleID between colData and proData is wrong please check your data")
  }

  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    return(
      data.frame(
        rowID = rownames(cormat)[row(cormat)[ut]],
        columnID = rownames(cormat)[col(cormat)[ut]],
        cor = (cormat)[ut],
        p = pmat[ut]
      )
    )
  }

  if(is.element(Object, "Feature")){
    # Features -> colnames; SampleID -> rownames
    cor <- Hmisc::rcorr(t(proData), type = "spearman")
  }else if(is.element(Object, "Meta")){
    colData_cln <- colData[, sapply(colData, function(x){is.numeric(x)})] %>% stats::na.omit()
    cor <- Hmisc::rcorr(as.matrix(colData_cln), type = "spearman")
  }else if(is.element(Object, "Mix")){
    mdat <- t(proData) %>% data.frame() %>%
      tibble::rownames_to_column("SampleID") %>%
      dplyr::inner_join(colData %>% tibble::rownames_to_column("SampleID"),
                        by = "SampleID")
    mdat_cln <- mdat[, sapply(mdat, function(x){is.numeric(x)})] %>% stats::na.omit()
    cor <- Hmisc::rcorr(as.matrix(mdat_cln), type = "spearman")
  }

  FCM <- flattenCorrMatrix(signif(cor$r, 5), signif(cor$P, 5))

  res <- list(corfit=cor, CorMatrix=FCM)
  return(res)
}
