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
#' library(convert)
#' data(ExprSet_species)
#' #data(ExprSet_species_count)
#' #data(ExprSet_genus)
#' #data(ExprSet_genus_count)
#' exprs(ExprSet_species)
#' pdata(ExprSet_species)
#' fData(ExprSet_species)
#'

# roxygen2::roxygenise("D:/Project/MyRtools/")
# roxygen2::roxygenise("~/Documents/Mygithub/MyRtools/")
# devtools::check(document = FALSE)
# devtools::build()
##################################################################
# building dataset
# library(dplyr)
# Profile <- data.table::fread(system.file("extdata", "Species_relative_abundance.tsv", package="MyRtools"))  %>% tibble::column_to_rownames("V1")
# Metadata <- read.csv(system.file("extdata", "Metadata.csv", package="MyRtools"))
# Feature <- read.csv(system.file("extdata", "Species_feature.csv", package="MyRtools")) %>% tibble::column_to_rownames("Species")
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
# save(ExprSetRawRB, file = "ExprSetRawRB.rda")
# Profile <- Profile * 10e6
# ExprSetRawCount <- get_ExprSet(profile=Profile,
#                                metadata=Metadata,
#                                feature=Feature,)
# save(ExprSetRawCount, file = "data/ExprSetRawCount.rda")
###########################################################################
