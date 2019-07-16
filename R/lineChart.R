#'  linechart1
#'
#' @description
#' show the linechart for data.
#'
#' @import ggplot2
#'
#' @return linechart
#' @export
#' @examples
#' linechart1()
#'
linechart1 <- function(){

  df <- data.frame(year=rep(1990:2015, times = 2),
                   type=rep(c('A','B'),each = 26),
                   value=c(runif(26),runif(26, min = 1,max = 1.5)))
  ggplot(data = df, mapping = aes(x = year, y = value, linetype = type,
                                  colour = type, shape = type, fill = type))+
    geom_line() + geom_point()+
    scale_linetype_manual(values = c(1,2))+
    scale_color_manual(values = c('steelblue','darkred'))+
    scale_shape_manual(values = c(21,23))+
    scale_fill_manual(values = c('red','black'))+
    theme_classic()
}

