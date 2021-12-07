#' @title Transform the expression profile in `assayData` sample by sample
#'
#' @description  Transform the taxa abundances in `assayData` sample by sample, which means
#' the counts of each sample will be transformed individually.
#'
#' @references https://github.com/yiluheihei/microbiomeMarker.
#'
#' @param object, Object; a [`matrix`] or [`assayData-class`] or [`ExpressionSet-class`].
#' @param transform, Character; transformation to apply, the options inclulde:
#' * "identity", return the original data without any transformation.
#' * "log10", the transformation is `log10(object)`, and if the data contains
#'   zeros the transformation is `log10(1 + object)`.
#' * "log10p", the transformation is `log10(1 + object)`.
#'
#' @return A object matches the class of argument `object` with the transformed
#'   `assayData`.
#'
#' @export
#' @examples
#'
transform_profile <- function(object,
                              transform = c("identity", "log10", "log10p")){

  transform <- match.arg(transform, c("identity", "log10", "log10p"))
  if(inherits(object, "ExpressionSet")){
    prf <- as(exprs(object), "matrix")
  }else if(inherits(object, "environment")){
    prf <- as(object$exprs, "matrix")
  }else{
    prf <- object
  }

  if (transform == "identity") {
    abd <- prf
  } else if (transform == "log10") {
    abd <- transform_log10(prf)
  } else {
    abd <- transform_log10p(prf)
  }

  if(inherits(object, "ExpressionSet")){
    object <- get_TransformedExprSet(object, abd)
  }else if(inherits(object, "environment")){
    object <- assayDataNew(exprs=abd)
  }else{
    object <- abd
  }

  return(object)
}

# the data is transformed using log10(1 + x) if the data contains zeroes
transform_log10 <- function(x) {
  if (min(x) == 0) {
    warning("Profile table contains zeroes. Using log10(1 + x) instead.")
    x_norm <- log10(1 + x)
  } else {
    x_norm <- log10(x)
  }

  return(x_norm)
}

# the data is transformed using log10(1 + x)
transform_log10p <- function(x) {
  log10(1 + x)
}
