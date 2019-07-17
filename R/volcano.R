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
volcano1 <- function(){

  # read.table("inst/results.txt", header = T) -> volcano_data
  # use_data(volcano_data)

  volcano_data$color <- ifelse(volcano_data$padj<0.05 & abs(volcano_data$log2FoldChange)>= 1,
                               ifelse(volcano_data$log2FoldChange > 1,'red','blue'),'gray')
  color <- c(red = "red",gray = "gray",blue = "blue")

  return(
    ggplot(volcano_data, aes(log2FoldChange, -log10(padj), col = color)) +
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


