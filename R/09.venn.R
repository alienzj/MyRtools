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
