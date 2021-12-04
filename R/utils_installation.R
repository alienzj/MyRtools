#' @title Installation of R packages
#'
#' @description
#' Installing the R packages coming from CRAN and bioconductor.
#'
#' @details 12/4/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param Package, Character; (Required).
#' @param Type, Character; the kind of repository(default: "CRAN").
#'
#' @return
#'
#' @export
#'
#' @usage InstallPackage(Package=Packages, Type="CRAN")
#' @examples
#'
#' \donttest{
#'    Packages_list <- c("DESeq2", "GSVA")
#'    InstallPackage(Package=Packages_list, Type="bioconductor")
#' }
#'
InstallPackage <- function(Package=packages, Type=c("CRAN", "bioconductor")){

  # Install Package not yet installed
  installed_packages <- Package %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    if(is.element(Type, "CRAN")){
      lapply(packages[!installed_packages], install.packages)
    }else{
      if(!require(BiocManager)){install.packages("BiocManager")}
      lapply(packages[!installed_packages], BiocManager::install)
    }
  }
  # Packages Loading
  invisible(lapply(packages, library, character.only = TRUE))

  message("The installation is finished")
}
