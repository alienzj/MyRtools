#'  heatmap1
#'
#' @description
#' show the heatmap for data.
#'
#' @import ggplot2
#' @importFrom tidyr gather
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
