#' @title Trim samples or features in `assayData` by Occurrences
#'
#' @description
#' trim samples or features in `assayData` by Occurrences,
#' which means the samples or features will be discard if they could not pass the cutoff.
#'
#' @details 12/6/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param object, Object; a [`matrix`] or [`assayData-class`] or [`ExpressionSet-class`].
#' @param cutoff, Numeric; the threshold for triming (default: Cutoff=0.2).
#' @param trim, Character; transformation to apply, the options inclulde:
#' * "none", return the original data without any trim.
#' * "both", trim the features and samples of profile more than cutoff.
#' * "feature", trim the features of profile more than cutoff.
#' * "sample", trim the samples of profile more than cutoff.
#' * "Group", trim the features per Group, the samples of the profile more than cutoff, respectively.
#'
#' @return
#'  A trimed `object` with profile more than cutoff.
#'
#' @export
#'
#' @importFrom Biobase exprs assayDataNew
#'
#' @usage run_trim(object, object=0.2, trim="Both")
#'
#' @examples
#' \donttest{
#'    data("ExprSet_species")
#'    run_trim(ExprSet_species, object=0.2, trim="both")
#' }
#'
run_trim <- function(object,
                     cutoff = 0.2,
                     trim = c("none", "both", "feature", "sample", "Group")){

  trim <- match.arg(trim, c("none", "both", "feature", "sample", "Group"))
  if(inherits(object, "ExpressionSet")){
    prf <- as(Biobase::exprs(object), "matrix")
  }else if(inherits(object, "environment")){
    prf <- as(object$exprs, "matrix")
  }else{
    prf <- object
  }

  if(trim == "feature"){
    tmp1 <- trim_FeatureOrSample(prf, 1, cutoff)
    remain_features <- rownames(tmp1)
    remain_samples <- colnames(prf)
  }else if(trim == "sample"){
    tmp2 <- trim_FeatureOrSample(prf, 2, cutoff)
    remain_features <- rownames(prf)
    remain_samples <- rownames(tmp2)
  }else if(trim == "both"){
    tmp1 <- trim_FeatureOrSample(prf, 1, cutoff)
    tmp2 <- trim_FeatureOrSample(prf, 2, cutoff)
    remain_features <- rownames(tmp1)
    remain_samples <- rownames(tmp2)
  }else if(all(trim == "Group", inherits(object, "ExpressionSet"))){
    tmp3 <- trim_eachGroup(object, trim, cutoff)
    remain_features <- rownames(tmp3$features)
    remain_samples <- rownames(tmp3$samples)
  }else if(trim == "none"){
    return(object)
  }
  prf_remain <- prf[remain_features, remain_samples]

  if(inherits(object, "ExpressionSet")){
    assayData(object) <- Biobase::assayDataNew(exprs=prf_remain)
    object <- get_TrimedExprSet(object)
  }else if(inherits(object, "environment")){
    object <- Biobase::assayDataNew(exprs=prf_remain)
  }else{
    object <- prf_remain
  }

  return(object)
}

# the data is trimmed by threshold
#' @keywords internal
trim_FeatureOrSample <- function(x, nRow, threshold){

  df_occ <- apply(x, nRow, function(x){
      length(x[c(which(!is.na(x)&x!=0))])/length(x)
    }) %>%
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
# the data is trimmed by threshold per group
#' @keywords internal
#' @importFrom dplyr %>% intersect select inner_join filter all_of
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom stats setNames
#' @importFrom Biobase exprs pData
#'
trim_eachGroup <- function(x, group_info, threshold){

  edata <- Biobase::exprs(x)
  pdata <- Biobase::pData(x)

  mdat <- pdata %>%  data.frame() %>%
    tibble::rownames_to_column("SampleID") %>%
    dplyr::select(dplyr::all_of(c("SampleID", group_info))) %>%
    dplyr::inner_join(edata %>% t() %>% data.frame() %>%
                        tibble::rownames_to_column("SampleID"),
            by = "SampleID")

  group_name <- unique(mdat$Group)
  group_name_each <- lapply(group_name, function(x){
     mdat_cln <- mdat %>% dplyr::filter(Group%in%x)
     return(mdat_cln$SampleID)
  })

  feature_occ_each <- sapply(1:length(group_name_each), function(i){
        df <- mdat %>% dplyr::filter(SampleID%in%group_name_each[[i]])
        df2 <- df %>% dplyr::select(-"Group") %>%
          tibble::column_to_rownames("SampleID") %>%t()
        ratios <- as.numeric(apply(df2, 1, function(x){length(x[c(which(!is.na(x)&x!=0))])/length(x)}))
        return(ratios)
  }) %>% data.frame() %>% stats::setNames(paste0(group_name, "_Occ"))

  rownames(feature_occ_each) <- colnames(mdat)[-c(1:2)]

  feature_KEEP <- apply(feature_occ_each > threshold, 1, all) %>%
    data.frame() %>% stats::setNames("Status") %>%
    dplyr::filter(Status)

  sample_occ <- apply(edata, 2, function(x){length(x[c(which(!is.na(x)&x!=0))])/length(x)}) %>%
    data.frame() %>% stats::setNames("Occ")
  sample_KEEP <- apply(sample_occ > threshold, 1, all) %>%
    data.frame() %>% stats::setNames("Status") %>%
    dplyr::filter(Status)

  res <- list(features=feature_KEEP,
              samples=sample_KEEP)

  return(res)

}
