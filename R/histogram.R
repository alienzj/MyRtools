#' histogram
#'
#' @description histogram visualization for daily work
#'
#' @import ggplot2
#' @return histogram
#'
#' @export
#'
#' @examples
#' histogram1()
#'
histogram1 <- function(){

  dat<-data.frame(x=rnorm(2000))

  return(
    ggplot(dat,aes(x,fill=cut(x,100)))+
      geom_histogram(bins=50,show.legend=FALSE)+
      scale_fill_discrete(h=c(250,10),c=120,l=70)+
      theme_minimal()+
      labs(x="VariableX",y="n")+
      ggtitle("HistogramofX",subtitle=R.version.string)+
      labs(caption="zsrnog")+
      theme(panel.background=element_rect(fill='black'),
            plot.background=element_rect(fill='black'),
            plot.title=element_text(colour="blue"),
            plot.subtitle=element_text(colour="blue"),
            plot.caption=element_text(colour="blue"),
            axis.line=element_line(colour="grey80"),
            axis.text=element_text(colour="blue"),
            axis.title=element_text(colour="grey80"))
  )
}
