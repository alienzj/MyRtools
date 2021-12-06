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
# devtools::check(document = FALSE)
# devtools::build()
##################################################################
# building dataset
# library(dplyr)
# Profile <- data.table::fread(system.file("extdata", "Species_relative_abundance.tsv", package="MyRtools"))  %>% tibble::column_to_rownames("V1")
# Metadata <- read.csv(system.file("extdata", "Metadata.csv", package="MyRtools"))
# Feature <- read.csv(system.file("extdata", "Species_feature.csv", package="MyRtools")) %>% tibble::column_to_rownames("Species")

# ExprSet_species <- get_ExprSet(profile=Profile,
#                                metadata=Metadata,
#                                feature=Feature,
#                                trim=TRUE,
#                                occ_Feature=0.2,
#                                each=TRUE,
#                                occ_Sample=0.2)
# save(ExprSet_species, file = "data/ExprSet_species.rda")
# Profile <- Profile * 10e6
# ExprSet_species_count <- get_ExprSet(profile=Profile,
#                                      metadata=Metadata,
#                                      feature=Feature,
#                                      trim=TRUE,
#                                      occ_Feature=0.2,
#                                      each=TRUE,
#                                      occ_Sample=0.2)
# save(ExprSet_species_count, file = "data/ExprSet_species_count.rda")
###########################################################################
# PROF <- Profile_table(as.matrix(Profile), feature_are_row=TRUE)
# META <- Sample_data(Metadata) %>% tibble::column_to_rownames("SampleID")
# FEATURE <- Feature_table(as.matrix(Feature))
#
# mydataset <- MyDataSet(Profile_table=PROF,
#                        Sample_data=META,
#                        Feature_table=FEATURE)
# save(mydataset, file = "data/MyDataSet.rda")
#
###########################################################################
# Profile <- data.table::fread(system.file("extdata", "Genus_relative_abundance.tsv", package="MyRtools"))  %>% tibble::column_to_rownames("V1")
# Metadata <- read.csv(system.file("extdata", "Metadata.csv", package="MyRtools"))
# Feature <- read.csv(system.file("extdata", "Genus_feature.csv", package="MyRtools")) %>% tibble::column_to_rownames("Genus")
#
#
# ExprSet_genus <- get_ExprSet(profile=Profile,
#                                metadata=Metadata,
#                                feature=Feature,
#                                trim=TRUE,
#                                occ_Feature=0.2,
#                                each=TRUE,
#                                occ_Sample=0.2)
# save(ExprSet_genus, file = "data/ExprSet_genus.rda")
# Profile <- Profile * 10e6
# ExprSet_genus_count <- get_ExprSet(profile=Profile,
#                                      metadata=Metadata,
#                                      feature=Feature,
#                                      trim=TRUE,
#                                      occ_Feature=0.2,
#                                      each=TRUE,
#                                      occ_Sample=0.2)
# save(ExprSet_genus_count, file = "data/ExprSet_genus_count.rda")
###########################################################################
