#' @title the new filtered ExpressionSet object
#'
#' @description
#' After filtering samples or features, we need to obtain the new ExpressionSet object.
#'
#' @details 12/2/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param dataset, ExpressionSet; (Required) unbalanced ExpressionSet object.
#'
#' @return
#' an filtered ExpressionSet object
#'
#' @export
#'
#' @importFrom dplyr %>% select all_of
#' @import Biobase
#'
#' @usage get_FilterExprSet(dataset=ExpressionSet)
#' @examples
#'
#' \donttest{
#'
#' }
#'
get_FilterExprSet <- function(dataset){

  FeatureNames <- rownames(exprs(dataset))
  SampleNames <- colnames(exprs(dataset))
  if(any(length(FeatureNames) != length(featureNames(dataset)),
         length(SampleNames) != length(sampleNames(dataset)))){

    pdata <- pData(dataset) %>% t() %>% data.frame() %>%
      dplyr::select(dplyr::all_of(SampleNames)) %>%
      t() %>% data.frame()
    if(!any(rownames(pdata) == SampleNames)){
      stop("Please check the order of SampleID between phen and prof")
    }
    pdata <- new("AnnotatedDataFrame", data=pdata)

    fdata <- fData(dataset) %>% t() %>% data.frame() %>%
      dplyr::select(dplyr::all_of(FeatureNames)) %>%
      t() %>% data.frame()
    if(!any(rownames(fdata) == FeatureNames)){
      stop("Please check the order of SampleID between phen and prof")
    }
    fdata <- new("AnnotatedDataFrame", data=fdata)
    experimentdata <- experimentData(dataset)

    res <- new("ExpressionSet",
               exprs=exprs(dataset),
               phenoData=pdata,
               featureData=fdata,
               experimentData=experimentdata)
  }else{
    res <- dataset
  }
  return(res)
}
################################################################################
#' @title the new transformed ExpressionSet object
#'
#' @description
#' After tranformed the expression, we need to obtain the new transformed ExpressionSet object.
#'
#' @details 12/2/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param dataset, ExpressionSet; (Required) old ExpressionSet object.
#' @param expr, matrix; (Required) the transformed expression matrix.
#'
#' @return
#' an transformed ExpressionSet object
#'
#' @export
#'
#' @importFrom dplyr %>% select all_of
#' @import Biobase
#'
#' @usage get_TransformedExprSet(dataset=ExpressionSet, expr=tranformedMatrix)
#' @examples
#'
#' \donttest{
#'
#' }
#'
get_TransformedExprSet <- function(dataset, expr){

  FeatureNames <- rownames(expr)
  SampleNames <- colnames(expr)
  if(!all(exprs(dataset) == expr)){

    pdata <- pData(dataset)
    if(!any(rownames(pdata) == SampleNames)){
      stop("Please check the order of SampleID between phen and prof")
    }
    pdata <- pData(dataset) %>% t() %>% data.frame() %>%
      dplyr::select(dplyr::all_of(SampleNames)) %>%
      t() %>% data.frame()
    pdata <- new("AnnotatedDataFrame", data=pdata)

    fdata <- fData(dataset)
    if(!any(rownames(fdata) == FeatureNames)){
      stop("Please check the order of SampleID between phen and prof")
    }
    fdata <- fData(dataset) %>% t() %>% data.frame() %>%
      dplyr::select(dplyr::all_of(FeatureNames)) %>%
      t() %>% data.frame()
    fdata <- new("AnnotatedDataFrame", data=fdata)

    experimentdata <- experimentData(dataset)

    res <- new("ExpressionSet",
               exprs=expr,
               phenoData=pdata,
               featureData=fdata,
               experimentData=experimentdata)
  }else{
    res <- dataset
  }
  return(res)
}
################################################################################