#'  pie1
#'
#' @description
#' show the pie for data
#'
#' @return pie
#' @export
#' @examples
#' pie1()
#'
pie1 <- function(){

  data <- read.csv(system.file("extdata", "vegetables.csv", package="patternplot"))
  pattern.type <- c('hdashes', 'vdashes', 'bricks')
  pattern.color <- c('red3','green3', 'white' )
  background.color <- c('dodgerblue', 'lightpink', 'orange')

  return(
    patternplot::patternpie(group=data$group,
               pct=data$pct,
               label=data$label,
               pattern.type=pattern.type,
               pattern.color=pattern.color,
               background.color=background.color,
               frame.color='grey40',
               pixel=0.3,
               pattern.line.size=0.3,
               frame.size=1.5,
               label.size=5,
               label.distance=1.35)+
      ggtitle('(B) Colors with Patterns')
  )
}
