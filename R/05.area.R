#'  area1
#'
#' @description
#' show the areaplot for data.
#'
#' @return areaplot
#' @export
#' @examples
#' area1()
#'
area1 <- function(dat=DATA){

  AREA <- dat %>% group_by(Fruit, Store) %>%
    summarize(mean_Weight=mean(Weight)) %>%
    ggplot(aes(x=Store, group=Fruit))+
      geom_area(aes(y=mean_Weight, fill=as.factor(Fruit)),
                position="stack",linetype="dashed")+
      geom_hline(aes(yintercept=mean(mean_Weight)),
                 color="blue",
                 linetype="dashed",
                 size=1)+
      guides(fill=guide_legend(title=NULL))+
      theme_bw(base_size=12) +
      theme(
        plot.title = element_text(size=10, color="black", face="bold", hjust=.5),
        axis.title = element_text(size=10, color="black", face="bold"),
        axis.text = element_text(size=9, color="black"),
        axis.ticks.length = unit(-0.05, "in"),
        axis.text.y = element_text(margin=unit(c(0.3,0.3,0.3,0.3), "cm"), size=9),
        axis.text.x = element_text(margin=unit(c(0.3,0.3,0.3,0.3), "cm")),
        text = element_text(size=8, color="black"),
        strip.text = element_text(size=9, color="black", face="bold"),
        panel.grid = element_blank())

  return(AREA)
}
