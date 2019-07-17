#'  corrplot1
#'
#' @description
#' show the corrplot for data.
#'
#' @return corrplot
#' @export
#' @examples
#' corrplot1()
#'
corrplot1 <- function(dat=DATA2){

  corrplot(corr=cor(dat[1:7]),
           order = "AOE",
           type="upper",
           tl.pos = "d")
  corrplot(corr = cor(dat[1:7]),
           add=TRUE,
           type="lower",
           method="number",
           order="AOE",
           diag=FALSE,
           tl.pos="n",
           cl.pos="n")
}
