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
#' @usage get_TrimedExprSet(dataset=ExpressionSet)
#' @examples
#'
#' \donttest{
#'
#' }
#'
get_TrimedExprSet <- function(dataset){

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
    experimentdata@other$notes <- "Trimed ExpressionSet"

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
#' After tranforming the expression, we need to obtain the new transformed ExpressionSet object.
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
    experimentdata@other$notes <- "Transformed ExpressionSet"
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
#' Returns the intersection of species and samples for the components of x
#'
#' This function is used internally as part of the infrastructure to ensure that
#' component data types in a ExpressionSet-object have exactly the same feature.
#' It relies heavily on the \code{\link{Reduce}} function to determine the
#' strictly common species.
#'
#' @usage intersect_feature(x)
#'
#' @param x (Required). A \code{\link{ExpressionSet-class}} object
#'  that contains 2 or more components
#'  that in-turn describe feature.
#'
#' @return Returns a character vector of only those feature that are present in
#'  all feature-describing components of \code{x}.
#'
#' @seealso \code{\link{Reduce}}, \code{\link{intersect}}
#' @keywords internal
#' @examples
#'
intersect_feature <- function(x){
  feature_vectors <- f_comp_es("featureNames", x)
  feature_vectors <- feature_vectors[!sapply(feature_vectors, is.null)]
  return( Reduce("intersect", feature_vectors) )
}
#' @keywords internal
intersect_samples <- function(x){
  sample_vectors <- f_comp_es("sampleNames", x)
  sample_vectors <- sample_vectors[!sapply(sample_vectors, is.null)]
  return( Reduce("intersect", sample_vectors) )
}
################################################################################
# A relatively fast way to access from ExpressionSet object components
# f - function name as character string
# expressionset - a ExpressionSet object (ExpressionSet-class instance)
#' @keywords internal
f_comp_es <- function(f, expressionset){
  sapply(names(getSlots("ExpressionSet")), function(i, ps){
    eval(parse(text=paste(f, "(es@", i, ")", sep="")))
  }, mydataset)
}
################################################################################
#' @title the new Normalized ExpressionSet object
#'
#' @description
#' After Normalizing the expression, we need to obtain the new Normalized ExpressionSet object.
#'
#' @details 12/2/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param dataset, ExpressionSet; (Required) old ExpressionSet object.
#' @param expr, matrix; (Required) the transformed expression matrix.
#'
#' @return
#' an Normalized ExpressionSet object
#'
#' @export
#'
#' @importFrom dplyr %>% select all_of
#' @import Biobase
#'
#' @usage get_NormalizedExprSet(dataset=ExpressionSet, expr=tranformedMatrix)
#' @examples
#'
#' \donttest{
#'
#' }
#'
get_NormalizedExprSet <- function(dataset, expr){

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
    experimentdata@other$notes <- "Normalized ExpressionSet"
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
#' @title the new imputed ExpressionSet object
#'
#' @description
#' After imputing features, we need to obtain the new ExpressionSet object.
#'
#' @details 12/13/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param dataset, ExpressionSet; (Required) imputed ExpressionSet object.
#'
#' @return
#' an imputed ExpressionSet object
#'
#' @export
#'
#' @importFrom dplyr %>% select all_of
#' @import Biobase
#'
#' @usage get_imputedExprSet(dataset=ExpressionSet)
#' @examples
#'
#' \donttest{
#'
#' }
#'
get_imputedExprSet <- function(dataset){

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
    experimentdata@other$notes <- "Imputed ExpressionSet"

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
#' @title Obtain the 95% CI Odd Ratio by glm model
#'
#' @details 12/13/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param datx, Data.frame; (Required) Metadata with SampleID and Group information.
#' @param daty, Data.frame; (Required) Profile with SampleID, ordered by datx (Row->Features; Column->SampleID).
#' @param GroupName, Character; (Required) the contrast group.
#'
#' @return
#' 95% CI Odd Ratio with featureID
#'
#' @export
#'
#' @importFrom dplyr %>% select all_of
#' @importFrom tibble rownames_to_column
#' @importFrom stats glm
#'
#' @usage run_OddRatio(datx, daty, GroupName)
#' @examples
#'
#' \donttest{
#'
#' }
#'
run_OddRatio <- function(datx, daty, GroupName){

  # glm result for odd ratios 95%CI
  mdat <- dplyr::inner_join(datx %>% tibble::rownames_to_column("SampleID") %>%
                              dplyr::select(dplyr::all_of(c("SampleID", "Group"))),
                            daty %>% t() %>% data.frame() %>%
                              tibble::rownames_to_column("SampleID"),
                            by = "SampleID") %>%
    tibble::column_to_rownames("SampleID")

  dat_phe <- mdat %>% dplyr::select(Group) %>%
    mutate(Group=ifelse(Group==GroupName[2], 1, 0))
  dat_prf <- mdat %>% dplyr::select(-Group)

  glmFun <- function(GroupN, MarkerN){

    MarkerN[MarkerN==0] <- min(MarkerN[MarkerN!=0])
    dat_glm <- data.frame(group=GroupN,
                          marker=scale(MarkerN, center=TRUE, scale=TRUE)) %>%
      na.omit()
    model <- summary(stats::glm(group ~ marker, data = dat_glm,
                                family = binomial(link = "logit")))
    res <- signif(exp(model$coefficients["marker", 1]) +
                    qnorm(c(0.025,0.5,0.975)) * model$coefficients["marker",1], 2)

    return(res)
  }

  glm_res <- t(apply(dat_prf, 2, function(x, group){
    res <- glmFun(group, as.numeric(x))
    return(res)
  }, dat_phe$Group))

  Odd <- glm_res %>% data.frame() %>%
    setNames(c("upper", "expected","lower")) %>%
    mutate("Odds Ratio (95% CI)" = paste0(expected, " (", lower, ";", upper, ")"))
  Odd$FeatureID <- rownames(glm_res)

  res <- Odd[, c(5, 4)]
  return(res)
}
