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
#' expres(ExprSet_species)
#' pdata(ExprSet_species)
#' fData(ExprSet_species)
#'


#' Example dataset
#'
#'
#' @title the Species integer data of ExpressionSet Object
#'
#' @description
#' A ExpressionSet Object from *get_ExprSet* function.
#'
#' @details 12/2/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @docType data
#' @name ExprSet_species_count
#' @format ExpressionSet Object contains Expression profile and Metadata information.
#' @source Generated from *get_ExprSet* function
#' @examples
#'
#' library(convert)
#' data(ExprSet_species_count)
#' expres(ExprSet_species_count)
#' pdata(ExprSet_species_count)
#' fData(ExprSet_species_count)
#'

# roxygen2::roxygenise("D:/Project/MyRtools/")
# devtools::check(document = FALSE)
# devtools::build()

# building dataset
# Profile <- data.table::fread(system.file("extdata", "Species_relative_abundance.tsv", package="MyRtools"))  %>% tibble::column_to_rownames("V1")
# Metadata <- read.csv(system.file("extdata", "Metadata.csv", package="MyRtools"))
# Feature <- read.csv(system.file("extdata", "Species_feature.csv", package="MyRtools")) %>% tibble::column_to_rownames("Species")
#
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
#
# PROF <- Profile_table(as.matrix(Profile), feature_are_row=TRUE)
# META <- Sample_data(Metadata)
# FEATURE <- Feature_table(Feature)
#
# mydataset <- MyDataSet(Profile_table=PROF,
#                        Sample_data=META,
#                        Feature_table=FEATURE)
# save(mydataset, file = "data/MyDataSet.rda")
#
