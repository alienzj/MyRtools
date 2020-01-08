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

#' Install suggest packages
#'
#' @description The install_suggest_packages function will install the suggest packages
#'
#' @param packages   a list of packages
#' @examples install_suggest_packages(packages_name_list)
#'

install_suggest_packages <- function(packages_name_list = c("randomForest")){
  usePackage <- function(p) {
    if (!is.element(p, installed.packages()[, 1]))
      install.packages(p, dep=TRUE, repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
    suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
  }
  invisible(lapply(packages_name_list, usePackage))
}


