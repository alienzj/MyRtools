#' Example dataset
#'
#'
#' @title the Species relative abundance of ExpressionSet Object
#'
#' @description
#' A ExpressionSet Object from *get_ExprSet* function.
#'
#' @details 12/2/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @docType data
#' @name ExprSet_species
#' @format ExpressionSet Object contains Expression profile and Metadata information.
#' @source Generated from *get_ExprSet* function
#' @examples
#'
#' library(Biobase)
#' data(ExprSetRawRB)
#' exprs(ExprSetRawRB)
#' pdata(ExprSetRawRB)
#' fData(ExprSetRawRB)
#'


################################
# Checking package
################################
# roxygen2::roxygenise("D:/Project/MyRtools/")
# roxygen2::roxygenise("~/Documents/Mygithub/MyRtools/")
# roxygen2::roxygenise("/disk/user/zouhua/github/MyRtools/")
# devtools::check(document = FALSE)
# devtools::build()


################################
# Building Dataset
################################
## loading data
# library(dplyr)
# Profile <- data.table::fread(system.file("extdata", "Species_relative_abundance.tsv", package="MyRtools"))  %>% tibble::column_to_rownames("V1")
# Metadata <- read.csv(system.file("extdata", "Metadata.csv", package="MyRtools"))
# Feature <- read.csv(system.file("extdata", "Species_feature.csv", package="MyRtools")) %>% tibble::column_to_rownames("Species")

## ExprSetRawRB
# Metadata$PID <- rep(paste0("P", formatC(1:100, digits = 1, flag = 0)), 3)
# Metadata <- Metadata %>% dplyr::select(PID, SampleID, everything())
# ExprSetRawRB <- get_ExprSet(profile=Profile,
#                        metadata=Metadata,
#                        feature=Feature,
#                        name="Hua Zou",
#                        lab="UCAS",
#                        contact="zouhua1@outlook.com",
#                        title="Experiment",
#                        abstract="Profile",
#                        url="www.zouhua.top",
#                        notes="Expression")
# save(ExprSetRawRB, file = "data/ExprSetRawRB.rda")

## ExprSetRawCount
# Profile <- Profile * 10e6
# ExprSetRawCount <- get_ExprSet(profile=Profile,
#                                metadata=Metadata,
#                                feature=Feature,)
# save(ExprSetRawCount, file = "data/ExprSetRawCount.rda")
