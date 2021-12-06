#' @title Transform the profile abundances in `Profile_table` sample by sample
#'
#' @description  Transform the taxa abundances in `Profile_table` sample by sample, which means
#' the counts of each sample will be transformed individually.
#'
#' @references https://github.com/yiluheihei/microbiomeMarker.
#'
#' @param object, Class; [`Profile_table-class`] or [`MyDataSet-class`].
#' @param transform, Character; transformation to apply, the options inclulde:
#' * "identity", return the original data without any transformation.
#' * "log10", the transformation is `log10(object)`, and if the data contains
#'   zeros the transformation is `log10(1 + object)`.
#' * "log10p", the transformation is `log10(1 + object)`.
#'
#' @return A object matches the class of argument `object` with the transformed
#'   `Profile_table`.
#'
#' @export
#' @examples
#' \donttest{
#' data(mydataset)
#' x1 <- transform_profile(mydataset)
#' head(Profile_table(x1), 10)
#' x2 <- transform_profile(mydataset, "log10")
#' head(Profile_table(x2), 10)
#' x3 <- transform_profile(mydataset, "log10p")
#' head(Profile_table(x3), 10)
#' }
#'
transform_profile <- function(object,
                                 transform = c("identity", "log10", "log10p")){

  transform <- match.arg(transform, c("identity", "log10", "log10p"))
  prf <- as(Profile_table(object), "matrix")

  if (transform == "identity") {
    abd <- prf
  } else if (transform == "log10") {
    abd <- transform_log10(prf)
  } else {
    abd <- transform_log10p(prf)
  }

  Profile_table(object) <- Profile_table(abd, feature_are_rows = feature_are_rows(object))

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
