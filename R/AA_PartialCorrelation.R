#' @title Association Analysis by Partial Correlation Coefficient
#'
#' @description
#'  The correlation between two variables (X, Y) while controling a third variable(Z).
#'  Here, we regard each Feature as X and each Metabolite as Y, the Age/Gender variable as Z.
#'
#' @details 12/4/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param Expression, ExpressionSet; (Required) ExpressionSet object.
#' @param Group_info, Character; design factor(default: "Group").
#' @param Group_name, Character; the group for filtering(default: NULL).
#' @param AdjVar, Character; (Required) the variables for correlation(default: AdjVar=c("Age", "Gender")).
#' @param Method, Character; (Required) the method for correlation(default: Method="pcor").
#'
#' @return
#'  a list of results
#'   Partial Correlation matrix(estimate, pvalue, statistics)
#'   Partial Correlation Coefficient table
#'
#' @export
#'
#' @importFrom dplyr %>% select filter intersect inner_join all_of
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom ppcor pcor.test spcor
#' @importFrom Biobase pData exprs
#'
#' @usage AA_PartialCC(dataset=ExpressionSet, Group_info="Group", Group_name=NULL, AdjVar=c("Age", "Gender"), Method="pcor")
#' @examples
#'
#' \donttest{
#' data(ExprSet_species)
#'
#' PartialCC_res <- AA_PartialCC(dataset=ExprSet_species, Group_info="Group", Group_name=NULL, AdjVar=c("Age", "Gender"), Method="pcor")
#'
#'}
#'
AA_PartialCC <- function(dataset=ExprSet_species,
                         Group_info="Group",
                         Group_name=NULL,
                         AdjVar=c("Age", "Gender"),
                         Method="pcor"){

  metadata <- Biobase::pData(dataset)
  colnames(metadata)[which(colnames(metadata) == Group_info)] <- "Group"
  profile <- Biobase::exprs(dataset)

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
    as.matrix()

  colVariable <- colData %>% dplyr::select(-dplyr::all_of(AdjVar))
  colVariable <- colVariable[, sapply(colVariable, function(x){is.numeric(x)})]

  Adjvariable <- colData %>% dplyr::select(dplyr::all_of(AdjVar))
  Adjvariable <- sapply(Adjvariable, function(x){
      if(is.numeric(x)){
        return(x)
      }else{
        return(as.numeric(factor(x)))
      }
    }
  )
  rownames(Adjvariable) <- rownames(colVariable)

  if(!all(rownames(colVariable) == colnames(proData))){
    stop("Order of sampleID between colVariable and proData is wrong please check your data")
  }
  num_row <- nrow(proData)
  num_col <- ncol(colVariable)
  matrix_estimate <- matrix(NA, nrow = num_row, ncol = num_col,
                            dimnames = list(rownames(proData),
                                            colnames(colVariable)))
  matrix_pvalue <- matrix(NA, nrow = num_row, ncol = num_col,
                          dimnames = list(rownames(proData),
                                          colnames(colVariable)))
  matrix_statistic <- matrix(NA, nrow = num_row, ncol = num_col,
                             dimnames = list(rownames(proData),
                                             colnames(colVariable)))
  cor_table <- data.frame()

  for(i in 1:num_row){
    for(j in 1:num_col){
      index <- which(!is.na(colVariable[, j]))
      Yvar <- colVariable[, j][index]
      Xvar <- proData[i, ][index]
      Zvar <- Adjvariable[index, ]
      if(is.element(Method, "pcor")){
        cor_res <- ppcor::pcor.test(Xvar, Yvar, Zvar, method = "spearman")
      }else if(is.element(Method, "spcor")){
        cor_res <- ppcor::spcor.test(Xvar, Yvar, Zvar, method = "spearman")
      }
      matrix_estimate[i, j] <- cor_res$estimate
      matrix_pvalue[i, j] <- cor_res$p.value
      matrix_statistic[i, j] <- cor_res$statistic
      cor_table <- rbind(data.frame(Xvar=rownames(proData)[i],
                                    Yvar=colnames(colVariable)[j],
                                    Zvar=paste(AdjVar, collapse = "_"),
                                    cor_res), cor_table)
    }
  }

  res <- list(cortable=cor_table, estimateMatrix=matrix_estimate,
              pvalueMatrix=matrix_pvalue, statisticMatrix=matrix_statistic)
  return(res)
}
