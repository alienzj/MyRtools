#' @title impute the expression profile in `assayData` by features
#'
#' @description
#' The missing values in expression profile maybe affect the statistical results, imputing the NAs or Zero values should taken into account.
#'
#' @details 12/13/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param object, Object; a [`matrix`] or [`assayData-class`] or [`ExpressionSet-class`].
#' @param normalize, Character; normalization to apply, the options inclulde:
#' * "none": return the original data without any transformation.
#' * "deletion": the deletion method is used when the probability of missing variable is same for all observations.
#' * "GlobalStandard": (Mean/ Mode/ Median Imputation), consists of replacing the missing data for a given attribute by the mean or median (quantitative attribute) or mode (qualitative attribute) of all known values of that variable.
#' * "KNN": the missing values of an attribute are imputed using the given number of attributes that are most similar to the attribute whose values are missing.
#'
#' @return
#' A object matches the class of argument `object` with the imputated profile.
#'
#' @export
#'
#' @importFrom Biobase exprs assayDataNew
#'
#' @usage run_impute(object, normalize=c("none", "deletion", "GlobalStandard", "KNN"))
#' @examples
#'
run_impute <- function(object,
                       impute = c("none", "deletion", "GlobalStandard", "KNN")){

  impute <- match.arg(impute, c("none", "GlobalStandard", "KNN"))
  if(inherits(object, "ExpressionSet")){
    prf <- as(Biobase::exprs(object), "matrix")
  }else if(inherits(object, "environment")){
    prf <- as(object$exprs, "matrix")
  }else{
    prf <- object
  }

  if(impute == "none"){
    abd <- prf
  }else if(impute == "deletion"){
    abd <- na.omit(prf)
  }else if(impute == "GlobalStandard"){
    abd <- impute_GlobalStandard(prf)
  }else if(impute == "KNN"){
    abd <- impute_knn(prf)
  }

  if(inherits(object, "ExpressionSet")){
    object <- get_TransformedExprSet(object, abd)
  }else if(inherits(object, "environment")){
    object <- Biobase::assayDataNew(exprs=abd)
  }else{
    object <- abd
  }

  return(object)
}
################################################################################
#' Global Standard imputation
#'
#' Global Standard imputation is used the Global Median and standard as the cutoff to impute the missing value.
#' We calculate the mean or median for all non missing values of that variable then replace missing value with mean or median
#'
#' @param object, Object; a [`matrix`](Row->Features; Column->Samples).
#'
#' @export
#' @rdname impute-methods
#' @aliases impute_GlobalStandard
#'
impute_GlobalStandard <- function(object,
                                  width=0.3,
                                  downshift=1.8){

  # Row->Features; Column->Samples
  if(length(object[is.na(object)]) > 1){
    set.seed(123)
    for(i in 1:ncol(object)){
      temp <- object[, i]
      temp_sd <- width * sd(temp, na.rm = TRUE)
      temp_mean <- mean(temp, na.rm = TRUE) - downshift * sd(temp, na.rm = TRUE)
      n_missing <- sum(is.na(temp))
      object[is.na(temp), i] <- rnorm(n_missing, mean = temp_mean, sd = temp_sd)
      object_imputed <- object
    }
  }else{
      object_imputed <- object
      message("No NAs were found\n")
  }
  return(object_imputed)
}
################################################################################
#' k-nearest neighborhood samples to impute the missing values
#'
#' @param object, Object; a [`matrix`](Row->Features; Column->Samples).
#'
#' @importFrom impute impute.knn
#'
#' @export
#' @rdname impute-methods
#' @aliases impute_knn
#'
impute_knn <- function(object){

  # Row->Features; Column->Samples
  if(length(object[is.na(object)]) > 1){
    fit <- impute::impute.knn(object, k=10)
    object_imputed <- fit$data
  }else{
    object_imputed <- object
    message("No NAs were found\n")
  }
  return(object_imputed)
}
################################################################################
