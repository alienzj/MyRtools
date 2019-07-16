#'  corrplot1
#'
#' @description
#' show the corrplot for data.
#'
#' @import corrplot
#'
#' @return corrplot
#' @export
#' @examples
#' corrplot1()
#'
corrplot1 <- function(){

  corrplot(corr=cor(mtcars[1:7]),
           order = "AOE",
           type="upper",
           tl.pos = "d")
  corrplot(corr = cor(mtcars[1:7]),
           add=TRUE,
           type="lower",
           method="number",
           order="AOE",
           diag=FALSE,
           tl.pos="n",
           cl.pos="n")
}
