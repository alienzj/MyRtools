#' boxplot1
#'
#' @description boxplot
#'
#' @param DATA Table.
#' @return boxplot
#' @export
#' @examples boxplot1(DATA)
#'
boxplot1 <- function(dat=DATA){

  pr <- unique(dat$Fruit)
  grp.col <- c("#999999", "#E69F00", "#56B4E9")
  BOXPLOT <- dat %>% mutate(Fruit=factor(Fruit)) %>%
    ggplot(aes(x=Fruit, y=Weight, color=Fruit))+
      stat_boxplot( geom='errorbar',
                    width=0.15)+
      geom_boxplot(aes(fill=Fruit), width=.4,
                   outlier.colour="black",
                   outlier.shape=21,
                   outlier.size=1)+
      stat_summary(fun.y=mean,
                   geom="point",
                   shape=16,
                   size=2,
                   color="black")+
      stat_summary(fun.data=function(x){
          return(data.frame(y=0.98*120,label=length(x)))
            },geom="text", hjust=0.5, color="red", size=6)+
      stat_compare_means(comparisons =
           list(c(pr[1], pr[2]),c(pr[1], pr[3]), c(pr[2], pr[3])),
                       label="p.signif",
                       method="wilcox.test")+
      labs(title="Weight of Fruit", x="Fruit", y="Weight (kg)")+
      scale_color_manual(values=grp.col,
                         labels=pr)+
      scale_fill_manual(values=grp.col,
                        labels=pr)+
      guides(color=F, fil=F)+
      scale_y_continuous(sec.axis = dup_axis(label = NULL,  name = NULL),
                         breaks = seq(90, 108, 2),
                         limits = c(90, 120))+
      theme_box()
  return(BOXPLOT)
}


#' boxplot2
#'
#' @description boxplot
#'
#' @param DATA Table.
#' @return boxplot
#' @export
#' @examples boxplot2(DATA)
#'
boxplot2 <- function(dat=DATA){

  BOXPLOT <- dat %>% mutate(Fruit=factor(Fruit)) %>%
    patternplot::patternboxplot(dat$Store, dat$Weight,
                   group=dat$Fruit,
                   pattern.type=c('nwlines', 'blank', 'waves'),
                   pattern.color=c('black','white', 'grey20'),
                   background.color=c('gold','lightpink', 'lightgreen'),
                   density=c(2, 1, 3),
                   pixel=1.2,
                   xlab='',
                   ylab='Weights, pounds',
                   legend.h=1.5,
                   legend.y.pos=0.5,
                   legend.ratio1=0.005,
                   legend.x.pos=0.15)+
      ggtitle('(B) Colors with Patterns')+
      theme_box()
  return(BOXPLOT)
}
