#' bar1
#'
#' @description a function for bar plot using ggplot2
#'
#' @import ggplot2
#' @import ggthemes
#' @importFrom tidyr gather
#'
#' @return bar1
#'
#' @export
#'
#' @examples
#' bar1()

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
#' @import patternplot
#'
#' @export
#'
#' @examples
#' bar2()


bar2 <- function() {

    data <- read.csv(system.file("extdata", "monthlyexp.csv", package = "patternplot"))
    data <- data[which(data$City == "City 1"), ]
    
    x <- factor(data$Type, c("Housing", "Food", "Childcare"))
    y <- data$Monthly_Expenses
    pattern.type <- c("hdashes", "blank", "crosshatch")
    pattern.color <- c("black", "black", "black")
    background.color <- c("white", "white", "white")
    density <- c(20, 20, 10)
    
    return(patternbar(data, x, y, group = NULL,
                      ylab = "Monthly Expenses, Dollar",
                      pattern.type = pattern.type,
                      pattern.color = pattern.color, 
                      background.color = background.color,
                      pattern.line.size = 0.5,
                      frame.color = c("black", "black", "black"),
                      density = density) + 
        ggtitle("(A) Black and White with Patterns"))
}


