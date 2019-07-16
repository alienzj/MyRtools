#' boxplot1
#'
#' @description boxplot by ggplot2
#'
#' @import ggplot2
#'
#' @return boxplot
#'
#' @export
#'
#' @examples boxplot1()
#'
boxplot1 <- function(){

  ToothGrowth$dose <- as.factor(ToothGrowth$dose)

  return(
    ggplot(ToothGrowth, aes(x=dose, y=len, color=dose))+
      geom_boxplot(outlier.colour="red", outlier.shape=7, outlier.size=4)+
      geom_dotplot(binaxis='y', stackdir='center', stackratio=1, dotsize=1)+
      stat_summary(fun.y=mean, geom="point", shape=23, size=4, color="black")+
      labs(title="Plotoflengthperdose",x="Dose(mg)",y="Length")+
      scale_color_manual(values=c("#999999","#E69F00","#56B4E9"))+
      scale_y_continuous(sec.axis = dup_axis(label = NULL,  name = NULL))+
      theme_bw()+
      theme(legend.position="right")
  )
}


#' boxplot2
#' @description boxplot
#'
#' @return boxplot
#'
#' @export
#'
#' @examples boxplot2()
#'
boxplot2 <- function(){

  d <- data.frame(riskScore = abs(rnorm(100)),
                  BMI = sample(1:2, 100, replace=T),
                  stage = sample(1:2, 100, replace=T),
                  age = sample(1:2, 100, replace=T),
                  gender = sample(1:2, 100, replace=T))

  myboxplot <- function(x, data, col = NULL, xlab, pvalue="auto") {
    boxplot(x, data, axes = FALSE, col = col)
    axis(1, at = 1:2, labels =FALSE)
    text(1:2, y=par()$usr[3]-0.08*(par()$usr[4]-par()$usr[3]),
         srt=60, xpd=T, adj=1, labels = xlab)
    if (pvalue == "auto") {
      pvalue <- round(t.test(x, data=data)$p.value, 3)
    }

    if (!is.null(pvalue)) {
      plab <- paste("p =", pvalue)
      text(1.5, y = par()$usr[4]*1.05, xpd=T, label=plab, col=col)
    }
  }

  layout(t(1:4))
  par(oma=c(2, 4, 4, 0), mar=c(5,2,1,1), cex=1)

  myboxplot(riskScore~age,
            data=d,
            col='red',
            xlab=c("age < 60", "age > 60"))
  axis(2, las=1)

  myboxplot(riskScore~gender,
            data=d,
            col='green',
            xlab=c("Male", "Female"))

  myboxplot(riskScore~stage,
            data=d,
            col='blue',
            xlab=c("pStage 1-2", "pStage 1-2"))

  myboxplot(riskScore~BMI,
            data=d,
            col='cyan',
            xlab=c("BMI < 24", "BMI > 24"))
}


#' boxplot3
#' @description boxplot
#'
#' @import patternplot
#'
#' @return boxplot
#'
#' @export
#'
#' @examples boxplot3()
#'
#'

boxplot3 <- function(){

  data <- read.csv(system.file("extdata", "fruits.csv", package="patternplot"))
  x <- data$Store
  y <- data$Weight
  group <- data$Fruit
  pattern.type <- c('nwlines', 'blank', 'waves')
  background.color <- c('white','gray80', 'white')
  pattern.color <- c('black','black', 'black')
  density <- c(2, 1, 3)
  background.color <- c('gold','lightpink', 'lightgreen')
  pattern.color <- c('black','white', 'grey20')

  return(
    patternboxplot(data,x, y,
                   group=group,
                   pattern.type=pattern.type,
                   pattern.color=pattern.color,
                   background.color=background.color,
                   density=density,
                   pixel=1.2,
                   xlab='',
                   ylab='Weights, pounds',
                   legend.h=1.5,
                   legend.y.pos=0.5,
                   legend.ratio1=0.005,
                   legend.x.pos=0.15)+
      ggtitle('(B) Colors with Patterns')
  )
}
