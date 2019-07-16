#'  area1
#'
#' @description
#' show the areaplot for data.
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom ggthemes theme_wsj
#'
#' @return areaplot
#' @export
#' @examples
#' area1()
#'
area1 <- function(){

  dat <- mpg %>% group_by(year,class) %>%
    summarize(mean_cty=mean(cty))

  ggplot(dat,aes(x=class,group=year))+
    geom_area(aes(y=mean_cty, fill=as.factor(year)),
              position="stack",linetype="dashed")+
    geom_hline(aes(yintercept=mean(mean_cty)),
               color="blue",
               linetype="dashed",
               size=1)+
    theme_wsj()+
    scale_fill_wsj()+
    guides(fill=guide_legend(title=NULL))
}
