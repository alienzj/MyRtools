#' @title Filter samples or features in `Profile_table` by Occurrences
#'
#' @description
#' Filter samples or features in `Profile_table` by Occurrences,
#' which means the samples or features will be discard if they could not pass the cutoff.
#'
#' @details 12/6/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param Object, Class; a [`Profile_table-class`] or [`MyDataSet-class`].
#' @param Cutoff, Numeric; the threshold for filtering (default: Cutoff=0.2).
#' @param Type, Character; the type of filtering data ("Feature" or "Sample").
#' @param Each, Logical; filtering Features each group or whole data(default: Each=FALSE).
#'
#' @return
#'  A filtered `object` with the Cutoff.
#'
#' @export
#'
#' @importFrom stats mad median quantile sd
#'
#' @usage NormalizeFun(Vector, type="Zscore")
#' @examples
#'
filter_profile <- function(object,
                           Cutoff = 0.2,
                           Type = c("Feature", "Sample", "Both")){
  # object = mydataset
  # Cutoff = 0.2
  # Type = "Feature"

  Type <- match.arg(Type, c("Feature", "Sample"))
  prf <- as(Profile_table(object), "matrix")

  if(Type == "Feature"){
    tmp <- trimFun(prf, 1, Cutoff)
    abd <- prf[rownames(prf)%in%rownames(tmp), ]
  }else if(Type == "Sample"){
    tmp <- trimFun(prf, 2, Cutoff)
    abd <- prf[, colnames(prf)%in%rownames(tmp)]
  }else if(Type == "Both"){
    tmp1 <- trimFun(prf, 1, Cutoff)
    tmp2 <- trimFun(prf, 2, Cutoff)
    abd <- prf[rownames(prf)%in%rownames(tmp),
               colnames(prf)%in%rownames(tmp)]
  }

  Profile_table(object) <- Profile_table(abd, feature_are_rows = feature_are_rows(object))

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
