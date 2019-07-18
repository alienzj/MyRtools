#' CorrMatrix
#'
#' @description The flattenCorrMatrix function will format the correlation matrix into a table of 4 columns:
#' row names, column names, the correlation coefficient between each variable and the others, and the p-values
#' cormat : matrix of the correlation coefficients
#' pmat : matrix of the correlation p-values
#'
#' @param DATA Table mtcars
#' @return Matrix
#' @export
#' @examples Correlation_matrix(DATA2)
#'

Correlation_matrix <- function(dat=DATA2){

  res <- Hmisc::rcorr(as.matrix(dat))

  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
      rowID = rownames(cormat)[row(cormat)[ut]],
      columnID = rownames(cormat)[col(cormat)[ut]],
      cor = (cormat)[ut],
      p = pmat[ut]
    )
  }

  return(flattenCorrMatrix(signif(res$r, 5), signif(res$P, 5)))
}


