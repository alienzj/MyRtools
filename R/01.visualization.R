#' scatterplot1
#'
#' @description scatterplot1
#'
#' @param DATA Table.
#' @return scatterplot
#' @export
#' @examples scatterplot1(DATA2)
#'
scatterplot1 <- function(dat=DATA2){

  SCATTER <- dat %>% mutate(cyl=factor(cyl)) %>%
  ggplot(aes(x=wt, y=mpg, shape=cyl, color=cyl))+
      geom_point(size=3, alpha=0.4)+
      geom_smooth(method=lm,  linetype="dashed",
                color="darkred", fill="blue")+
      geom_text(aes(label=rownames(dat)), size=4)+
      theme_bw(base_size=12)+
      theme(
        plot.title = element_text(size=10, color="black", face="bold", hjust=.5),
        axis.title = element_text(size=10, color="black", face="bold"),
        axis.text = element_text(size=9, color="black"),
        axis.ticks.length = unit(-0.05, "in"),
        axis.text.y = element_text(margin=unit(c(0.3,0.3,0.3,0.3), "cm"), size=9),
        axis.text.x = element_blank(),
        text = element_text(size=8, color="black"),
        strip.text = element_text(size=9, color="black", face="bold"),
        panel.grid = element_blank())

    return(SCATTER)

}


#' histogram
#'
#' @description histogram visualization for daily work
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


#' bar1
#'
#' @description a function for bar plot using ggplot2
#'
#' @return bar1
#'
#' @export
#'
#' @examples
#' bar1()
#'
bar1 <- function() {

  data <- data.frame(Conpany = c("Apple", "Google", "Facebook", "Amozon", "Tencent"),
                     Sale2013 = c(5000, 3500, 2300, 2100, 3100),
                     Sale2014 = c(5050, 3800, 2900, 2500, 3300),
                     Sale2015 = c(5050, 3800, 2900, 2500, 3300),
                     Sale2016 = c(5050, 3800, 2900, 2500, 3300))
  mydata <- tidyr::gather(data, Year, Sale, -Conpany)

  return(
    ggplot(mydata, aes(Conpany, Sale, fill = Year)) +
      geom_bar(stat = "identity", position = "dodge") +
      guides(fill = guide_legend(title = NULL)) +
      ggtitle("The Financial Performance of Five Giant") +
      scale_fill_wsj("rgby", "") +
      theme_wsj() +
      theme(axis.ticks.length = unit(0.5, "cm"),
            axis.title = element_blank())
  )
}


#' bar2
#' @description bar
#'
#' @return bar
#'
#' @export
#'
#' @examples
#' bar2()
#'
bar2 <- function() {

  data <- read.csv(system.file("extdata", "monthlyexp.csv", package = "patternplot"))
  data <- data[which(data$City == "City 1"), ]

  x <- factor(data$Type, c("Housing", "Food", "Childcare"))
  y <- data$Monthly_Expenses
  pattern.type <- c("hdashes", "blank", "crosshatch")
  pattern.color <- c("black", "black", "black")
  background.color <- c("white", "white", "white")
  density <- c(20, 20, 10)

  return(
    patternplot::patternbar(data, x, y, group = NULL,
                            ylab = "Monthly Expenses, Dollar",
                            pattern.type = pattern.type,
                            pattern.color = pattern.color,
                            background.color = background.color,
                            pattern.line.size = 0.5,
                            frame.color = c("black", "black", "black"),
                            density = density) +
      ggtitle("(A) Black and White with Patterns"))
}


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

  return(BOXPLOT)
}


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


#'  heatmap1
#'
#' @description
#' show the heatmap for data.
#'
#' @return heatmap
#' @export
#' @examples
#' heat1()
#'
heat1 <- function(){

  data <- as.data.frame(matrix(rnorm(9*10),9,10))
  rownames(data) <- paste("Gene",1:9,sep="_")
  colnames(data) <- paste("sample",1:10,sep="_")
  data$ID <- rownames(data)
  data_m <- tidyr::gather(data, sampleID, value, -ID)

  ggplot(data_m,aes(x=sampleID,y=ID))+
    geom_tile(aes(fill=value))+
    scale_fill_gradient2("Expression",
                         low="green",
                         high="red",
                         mid="black")+
    xlab("samples")+
    theme_classic()+
    theme(axis.ticks=element_blank(),
          axis.line=element_blank(),
          panel.grid.major=element_blank(),
          legend.key=element_blank(),
          axis.text.x=element_text(angle=45,hjust=1,vjust=1),
          legend.position="top")
}


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


#'  venn1
#'
#' @description
#' show the venn for data.
#'
#' @return venn
#' @export
#' @examples
#' venn1()
#'
venn1 <- function() {

  A <- sample(LETTERS, 18, replace = FALSE)
  B <- sample(LETTERS, 18, replace = FALSE)
  C <- sample(LETTERS, 18, replace = FALSE)
  D <- sample(LETTERS, 18, replace = FALSE)

  venn.diagram(
    x = list(A = A, D = D, B = B, C = C),
    filename = "Group4.png",
    height = 450, width = 450,
    resolution = 300, imagetype = "png",
    col = "transparent",
    fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
    alpha = 0.5, cex = 0.45, cat.cex = 0.45)
}


#'  venn2
#'
#' @description
#' show the linechart for data.
#'
#' @return venn2
#' @export
#' @examples
#' venn2()
#'
venn2 <- function() {

  movies <- read.csv(system.file("extdata", "movies.csv",
                                 package = "UpSetR"), header = T, sep = ";")
  mutations <- read.csv(system.file("extdata", "mutations.csv",
                                    package = "UpSetR"), header = T, sep = ",")

  another.plot <- function(data, x, y) {
    round_any_new <- function(x, accuracy, f=round){f(x/ accuracy) * accuracy}
    data$decades <- round_any_new(as.integer(unlist(data[y])), 10, ceiling)
    data <- data[which(data$decades >= 1970), ]
    myplot <- (ggplot(data, aes_string(x = x)) +
                 geom_density(aes(fill = factor(decades)), alpha = 0.4) +
                 theme_bw()+
                 theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
                       legend.key.size = unit(0.4, "cm"))
    )
  }

  return(
    upset(movies,
          main.bar.color = "black",
          mb.ratio = c(0.5, 0.5),
          queries = list(
            list(query = intersects, params = list("Drama"), color = "red", active = F),
            list(query = intersects, params = list("Action", "Drama"), active = T),
            list(query = intersects, params = list("Drama", "Comedy", "Action"), color = "orange", active = T)),
          attribute.plots = list(gridrows = 50,
                                 plots = list(
                                   list(plot = histogram, x = "ReleaseDate", queries = F),
                                   list(plot = scatter_plot, x = "ReleaseDate", y = "AvgRating", queries = T),
                                   list(plot = another.plot, x = "AvgRating", y = "ReleaseDate", queries = F)
                                 ),
                                 ncols = 3)
    )
  )
}


#'  volcano1
#'
#' @description
#' show the volcano for data.
#'
#' @import ggplot2
#'
#' @return volcano
#' @export
#' @examples
#' volcano1()
#'
volcano1 <- function(dat=DATA3){

  dat$color <- with(dat, ifelse(padj<0.05 & abs(log2FoldChange)>= 1,
                                ifelse(log2FoldChange > 1,'red','blue'),'gray'))
  color <- c(red = "red",gray = "gray",blue = "blue")

  return(
    ggplot(dat, aes(log2FoldChange, -log10(padj), col = color)) +
      geom_point() +
      theme_bw() +
      scale_color_manual(values = color) +
      labs(x="log2 (fold change)",y="-log10 (q-value)") +
      geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
      geom_vline(xintercept = c(-1, 1), lty=4,col="grey",lwd=0.6) +
      theme(legend.position = "none",
            panel.grid=element_blank(),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14))
  )
}


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



#' circos1
#' @description circos
#'
#' @return circos
#'
#' @export
#'
#' @examples
#' circos1()
#'
circos1 <- function(){
  my.data<- matrix(
    c( 674,441,289,146,260,250,346,173,246,168,224,241,89 ,165,254,205,274,72,432,145, 70,119,189, 1,
       1307,893,793,487,442,468,652,392,555,467,663,620,220,413,487,465,657,91,856,264,124,274,359,10,
       388,364,154,144,160,134,228,100,211,115,195,152,75 ,196,122,124,204,49,178,133, 57,85 ,111, 3,
       389,312,177,182,171,145,303,94 ,272,147,204,158,74 ,134,110,181,207,68,237,128, 53,101,85 ,14,
       313,285,181,104,120,92 ,212,110,173,91 ,241,158,62 ,141,126,119,184,49,213,89 , 46,79 ,63 , 3,
       409,365,241,165,226,193,269,103,251,117,249,216,96 ,105,177,159,238,60,242,129, 61,99 ,62 , 8),6,24)
  rownames(my.data) <- c("A11","A12","A15","A16","A17","A18")
  colnames(my.data) <- c("1","2","3","4","5","6",
                         "7","8","9","10","11","12",
                         "13","14","15","16","17","18",
                         "19","20","21","22","X","Y")

  grid.col = NULL
  grid.col[c("A11", "A12", "A15", "A16","A17","A18")] = c("red", "yellow","green", "blue","purple","orange")
  grid.col[c("1","2","3","4","5","6",
             "7","8","9","10","11","12",
             "13","14","15","16","17","18",
             "19","20","21","22","X","Y")]="grey"

    circos.par(start.degree = 60,gap.degree = c(rep(2, nrow(my.data)-1), 8, rep(2, ncol(my.data)-1), 8))
    return(
      chordDiagram(my.data,
                   directional = FALSE,
                   diffHeight = 0.06,
                   grid.col = grid.col,
                   transparency = 0.5)
    )

}


