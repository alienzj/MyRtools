#' @title Filter samples or features in `Profile_table` by Occurrences
#'
#' @description
#' Filter samples or features in `Profile_table` by Occurrences,
#' which means the samples or features will be discard if they could not pass the cutoff.
#'
#' @details 12/6/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param Object, Class; a [`MyDataSet-class`].
#' @param Cutoff, Numeric; the threshold for filtering (default: Cutoff=0.2).
#' @param Type, Character; the type of filtering data ("Both", "Features", "Samples").
#'
#' @return
#'  A filtered `object` with the Cutoff.
#'
#' @export
#'
#' @importFrom stats mad median quantile sd
#'
#' @usage filter_method(Object, Cutoff=0.2, Type="Both")
#' @examples
#' data(mydataset)
#' filter_method(Object, Cutoff=0.2, Type="Both")
#'
filter_method <- function(object,
                          Cutoff = 0.2,
                          Type = c("Both", "Features", "Samples")){

  Type <- match.arg(Type, c("Both", "Features", "Samples"))
  prf <- as(Profile_table(object), "matrix")

  if(Type == "Features"){
    tmp1 <- trimFun(prf, 1, Cutoff)
    remain_features <- rownames(tmp1)
    remain_samples <- colnames(prf)
  }else if(Type == "Samples"){
    tmp2 <- trimFun(prf, 2, Cutoff)
    remain_features <- rownames(prf)
    remain_samples <- rownames(tmp2)
  }else if(Type == "Both"){
    tmp1 <- trimFun(prf, 1, Cutoff)
    tmp2 <- trimFun(prf, 2, Cutoff)
    remain_features <- rownames(tmp1)
    remain_samples <- rownames(tmp2)
  }
  object <- prune_feature(remain_features, object)

  # Verify there is more than one component
  # that describes samples before attempting to reconcile.
  object <- prune_samples(remain_samples, object)

  # Force both samples and feature indices to be in the same order.
  object <- index_reorder(object, "Both")

  #object@Profile_table <- Profile_table(abd, feature_are_rows = feature_are_rows(object))

  object
}

# the data is trimmed by threshold
trimFun <- function(x, nRow, threshold){

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
