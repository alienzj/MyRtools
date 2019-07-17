#' linechart1
#'
#' @description linechart
#'
#' @param DATA Table
#' @return linechart
#' @export
#' @examples linechart1(DATA)
#'
linechart1 <- function(dat=DATA){

  pr <- unique(dat$Fruit)
  grp.col <- c("#999999", "#E69F00", "#56B4E9")
  dat.cln <- sampling::strata(dat, stratanames="Fruit",
                   size = rep(round(nrow(dat)*0.1/3, -1), 3),
                   method="srswor")

  LINE <- dat %>% slice(dat.cln$ID_unit) %>%
    mutate(Year=as.character(rep(1996:2015, times = 3))) %>%
    mutate(Year=factor(as.character(Year))) %>%

    ggplot(aes(x=Year, y=Weight, linetype=Fruit, colour=Fruit,
               shape=Fruit, fill=Fruit))+
      geom_line(aes(group=Fruit)) +
      geom_point()+
      scale_linetype_manual(values = c(1:3))+
      scale_shape_manual(values = c(19, 21, 23))+
      scale_color_manual(values=grp.col,
                         labels=pr)+
      scale_fill_manual(values=grp.col,
                        labels=pr)+
      theme_bw()+
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

  return(LINE)
}
