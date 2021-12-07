#' @title Filter samples or features in `assayData` by Occurrences
#'
#' @description
#' Filter samples or features in `assayData` by Occurrences,
#' which means the samples or features will be discard if they could not pass the cutoff.
#'
#' @details 12/6/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param Object, Object; a [`matrix`] or [`assayData-class`] or [`ExpressionSet-class`].
#' @param Cutoff, Numeric; the threshold for filtering (default: Cutoff=0.2).
#' @param FilterType, Character; the type of filtering data ("Both", "Features", "Samples").
#'
#' @return
#'  A filtered `object` with the Cutoff.
#'
#' @export
#'
#' @importFrom stats mad median quantile sd
#' @import Biobase
#'
#' @usage filter_method(Object, Cutoff=0.2, FilterType="Both")
#'
#' @examples
#' \donttest{
#'    data("ExprSet_species")
#'    filter_method(ExprSet_species, Cutoff=0.2, FilterType="Both")
#' }
#'
filter_method <- function(object,
                          Cutoff = 0.2,
                          FilterType = c("Both", "Features", "Samples")){

  # object=ExprSet_species
  # Cutoff = 0.5
  # FilterType = "Both"

  FilterType <- match.arg(FilterType, c("Both", "Features", "Samples"))
  if(inherits(object, "ExpressionSet")){
    prf <- as(exprs(object), "matrix")
  }else if(inherits(object, "environment")){
    prf <- as(object$exprs, "matrix")
  }else{
    prf <- object
  }

  if(FilterType == "Features"){
    tmp1 <- trim_FeatureOrSample(prf, 1, Cutoff)
    remain_features <- rownames(tmp1)
    remain_samples <- colnames(prf)
  }else if(FilterType == "Samples"){
    tmp2 <- trim_FeatureOrSample(prf, 2, Cutoff)
    remain_features <- rownames(prf)
    remain_samples <- rownames(tmp2)
  }else if(FilterType == "Both"){
    tmp1 <- trim_FeatureOrSample(prf, 1, Cutoff)
    tmp2 <- trim_FeatureOrSample(prf, 2, Cutoff)
    remain_features <- rownames(tmp1)
    remain_samples <- rownames(tmp2)
  }
  prf_remain <- prf[remain_features, remain_samples]

  if(inherits(object, "ExpressionSet")){
    assayData(object) <- assayDataNew(exprs=prf_remain)
    object <- get_FilterExprSet(object)
  }else if(inherits(object, "environment")){
    object <- assayDataNew(exprs=prf_remain)
  }else{
    object <- prf_remain
  }

  return(object)
}

# the data is trimmed by threshold
#' @keywords internal
trim_FeatureOrSample <- function(x, nRow, threshold){

  df_occ <- apply(x, nRow, function(x){length(x[x!=0])/length(x)}) %>%
    data.frame() %>% stats::setNames("Occ") %>%
    tibble::rownames_to_column("type")
  if(nRow == 1){
    rownames(df_occ) <- rownames(x)
  }else{
    rownames(df_occ) <- colnames(x)
  }
  df_KEEP <- apply(df_occ > threshold, 1, all) %>%
    data.frame() %>% stats::setNames("Status") %>%
    dplyr::filter(Status)

  return(df_KEEP)
}
