#'  scatterplot1
#'
#' @description
#' show the scatterplot for data.
#'
#' @import gcookbook
#' @import ggplot2
#' @importFrom ggthemes theme_wsj
#'
#' @return scatterplot
#' @export
#' @examples
#' scatterplot1()
#'
scatterplot1 <- function(){

  return(
    ggplot(heightweight, aes(x=ageYear, y=heightIn,shape=sex))+
      geom_point(size=3, alpha=0.4)+
      stat_smooth(method=lm, level=0.99)+
      geom_text(aes(label=sex), size=4)+
      theme_wsj()
  )
}


