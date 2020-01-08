#' PERMANOVA
#'
#' @description PERMANOVA
#' PERMANOVA the over gut microbial species and each phenotype
#' @details 08/01/2020  ShenZhen China
#' @author  Hua Zou
#'
#' @usage PERMANOVA(phen, prof)
#' @examples  per_res <- PERMANOVA(phen, prof)
#'
PERMANOVA <- function(x, y) {

  sid <- intersect(as.character(x$SampleID), colnames(y))
  phe <- x %>% filter(SampleID %in% sid)
  prf <-  y %>% select(as.character(phe$SampleID)) %>%
    t() %>% data.frame()
  per <- apply(phe %>% select(-one_of("SampleID")), 2, function(a, pf){
    dat <- data.frame(value = a, pf)
    datphe <- dat$value %>% varhandle::unfactor()
    if (length(datphe) == 0 | unique(datphe) == 1) {
      res <- data.frame(length(datphe), rep(NA, 6))
      next
    }
    if (length(unique(datphe)) < 6) {
      datphe <- as.factor(datphe)
    }
    datprf <- dat[, -1, F]
    dis <- vegdist(datprf, method = "bray")
    set.seed(123)
    ad <- adonis(dis ~ datphe, permutations = 1000)
    tmp <- as.data.frame(ad$aov.tab) %>% slice(1)
    res <- c(length(datphe), as.numeric(tmp[, c(1:6)]))
    return(res)
  }, prf) %>% t() %>% data.frame()

  colnames(per) <- c("SumsOfSample", "Df", "SumsOfSqs",
                     "MeanSqs", "F.Model", "R2", "Pr(>F)")
  per$FDR <- p.adjust(per$`Pr(>F)`, method = "BH")
  return(per)
}


