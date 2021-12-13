#' @title Permutational Multivariate Analysis of Variance (PERMANOVA) Using Distance Matrices
#'
#' @description
#' The PERMANOVA is aimed to assess the association between the overall profile and metadata information.
#' It comprises permutation and ANOVA by using the F distribution to acquire the statistic.
#'
#' @details 12/4/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param Expression, ExpressionSet; (Required) A \code{ExpressionSet} object contains expression profile and
#' metadata including the continuous and categorical variables.
#' @param Group_info, Character; design factor(default: "Group").
#' @param Group_name, Character; the group for filtering(default: NULL).
#' @param Distance, Character; Normalizing feature(default: Distance="bray").
#'
#' @return
#'   PERMANOVA's results
#'
#' @export
#'
#' @importFrom dplyr %>% select all_of slice filter
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom vegan adonis vegdist
#' @importFrom stats setNames
#' @importFrom Biobase pData exprs
#'
#' @usage run_PERMANOVA(dataset=ExpressionSet, Group_info="Group", Group_name=NULL, Distance="bray")
#' @examples
#'
#' \donttest{
#' data("ExprSetRawRB")
#'
#' PERMANOVA_res <- run_PERMANOVA(dataset=ExprSetRawRB, Group_info="Group", Group_name=NULL, Distance="bray")
#' }
#'
run_PERMANOVA <- function(dataset=ExprSet,
                          Group_info="Group",
                          Group_name=NULL,
                          Distance="bray"){

  metadata <- Biobase::pData(dataset)
  colnames(metadata)[which(colnames(metadata) == Group_info)] <- "Group"
  profile <- Biobase::exprs(dataset)

  if(any(profile < 0) & is.element(Distance, "bray")){
    message("The value of profile is negative and the approach to calculate distance beteween samples shouldn't be bray crutis")
    Distance <- "euclidean"
  }

  # Prepare for input data
  if(is.null(Group_name)){
    colData <- metadata
  }else{
    phen <- metadata %>% dplyr::filter(Group%in%Group_name)
    intersect_sid <- dplyr::intersect(rownames(phen), colnames(profile))
    colData <- phen %>% tibble::rownames_to_column("SampleID") %>%
      dplyr::filter(SampleID%in%intersect_sid) %>%
      tibble::column_to_rownames("SampleID")
  }

  proData <- profile %>% data.frame() %>%
    dplyr::select(dplyr::all_of(rownames(colData))) %>%
    t()

  # distance
  sample_dis_all <- vegan::vegdist(proData, method = Distance)

  if(!all(rownames(colData) == rownames(proData))){
    stop("Order of sampleID between colData and proData is wrong please check your data")
  }


  perfun <- function(x, y, distance){
    index <- which(!is.na(x))
    datphe <- x[index]
    if(any(length(datphe) == 0, unique(datphe) == 1)){
      res <- data.frame(length(datphe), rep(NA, 6))
    }
    if(length(unique(datphe)) < 10){  # factor the variable if the levels of variable is less than 10
      datphe <- factor(datphe)
    }

    #set.seed(123)
    if(length(index) != length(datphe)){
      datprf <- y[index, ]
      sample_dis <- vegan::vegdist(datprf, method = Distance)
    }else{
      sample_dis <- distance
    }

    ad <- vegan::adonis(sample_dis ~ datphe, permutations = 999)
    tmp <- data.frame(ad$aov.tab) %>% dplyr::slice(1)
    res <- c(length(datphe), as.numeric(tmp[, c(1:6)]))

    return(res)
  }

  # numeric
  colData_numeric <- colData[, sapply(colData, function(x){is.numeric(x)})]
  # non-numeric
  colData_non_numeric <- colData[, sapply(colData, function(x){!is.numeric(x)})]

  if(ncol(colData_numeric) > 1){
    res_numeric  <- apply(colData_numeric, 2, perfun, proData, sample_dis_all) %>%
      t() %>% data.frame() %>%
      tibble::rownames_to_column("Phenotype")
  }else{
    res_numeric <- data.frame(rep(NA, 8))
  }

  if(ncol(colData_non_numeric) > 1){
    res_non_numeric  <- apply(colData_non_numeric %>% tibble::rownames_to_column("SampleID"),
                          2, perfun, proData, sample_dis_all) %>%
      t() %>% data.frame() %>%
      tibble::rownames_to_column("Phenotype")
  }else{
    res_non_numeric <- data.frame(rep(NA, 8))
  }

  res <- rbind(res_non_numeric, res_numeric)


  colnames(res) <- c("Phenotype", "SumsOfSample", "Df", "SumsOfSqs",
                     "MeanSqs", "F.Model", "R2", "Pr(>F)")

  res$FDR <- p.adjust(res$`Pr(>F)`, method = "BH")

  return(res)
}
