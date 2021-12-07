#' @title Filter samples or features in `assayData` by Occurrences
#'
#' @description
#' Filter samples or features in `assayData` by Occurrences,
#' which means the samples or features will be discard if they could not pass the cutoff.
#'
#' @details 12/6/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param object, Object; a [`matrix`] or [`assayData-class`] or [`ExpressionSet-class`].
#' @param cutoff, Numeric; the threshold for filtering (default: Cutoff=0.2).
#' @param filterType, Character; the type of filtering data ("identity", "both", "feature", "sample").
#'
#' @return
#'  A filtered `object` with the cutoff
#'
#' @export
#'
#' @importFrom stats mad median quantile sd
#' @import Biobase
#'
#' @usage filter_method(object, object=0.2, filterType="Both")
#'
#' @examples
#' \donttest{
#'    data("ExprSet_species")
#'    filter_method(ExprSet_species, object=0.2, filterType="both")
#' }
#'
filter_method <- function(object,
                          cutoff = 0.2,
                          filterType = c("identity", "both", "feature", "sample")){

  filterType <- match.arg(filterType, c("identity", "both", "feature", "sample"))
  if(inherits(object, "ExpressionSet")){
    prf <- as(exprs(object), "matrix")
  }else if(inherits(object, "environment")){
    prf <- as(object$exprs, "matrix")
  }else{
    prf <- object
  }

  if(filterType == "feature"){
    tmp1 <- trim_FeatureOrSample(prf, 1, cutoff)
    remain_features <- rownames(tmp1)
    remain_samples <- colnames(prf)
  }else if(filterType == "sample"){
    tmp2 <- trim_FeatureOrSample(prf, 2, cutoff)
    remain_features <- rownames(prf)
    remain_samples <- rownames(tmp2)
  }else if(filterType == "both"){
    tmp1 <- trim_FeatureOrSample(prf, 1, cutoff)
    tmp2 <- trim_FeatureOrSample(prf, 2, cutoff)
    remain_features <- rownames(tmp1)
    remain_samples <- rownames(tmp2)
  }else if(filterType == "identity"){
    return(object)
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
