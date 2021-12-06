################################################################################
# https://github.com/joey711/phyloseq/blob/master/R/allClasses.R
#' The S4 class for storing Profile_table information.
#'
#' Because orientation of these tables can vary by method, the orientation is
#' defined explicitly in the \code{feature_are_rows} slot (a logical).
#' The \code{Profile_table} class inherits the \code{\link{matrix}} class to store
#' expression values.
#' Various standard subset and assignment nomenclature has been extended to apply
#' to the \code{Profile_table} class, including square-bracket, \code{\link{t}}, etc.
#'
#' \describe{
#'    \item{feature_are_rows}{
#'		A single logical specifying the orientation of the expression table.
#'    }
#'
#'    \item{.Data}{This slot is inherited from the \code{\link{matrix}} class.}
#'  }
#' @name Profile_table-class
#' @rdname Profile_table-class
#' @exportClass Profile_table
setClass("Profile_table", representation(feature_are_rows="logical"), contains = "matrix")
################################################################################
#' The S4 for storing sample variables.
#'
#' Row indices represent samples, while column indices represent experimental
#' categories, variables (and so forth) that describe the samples.
#'
#' \describe{
#'
#'    \item{.Data}{data-frame data, inherited from the data.frame class.}
#'
#'    \item{row.names}{
#'	     Also inherited from the data.frame class;
#'       it should contain the sample names.
#'    }
#'
#'    \item{names}{Inherited from the data.frame class.}
#'
#'  }
#'
#' @name Sample_data-class
#' @rdname Sample_data-class
#' @exportClass Sample_data
setClass("Sample_data", contains="data.frame")
################################################################################
#' An S4 class that holds feature attributions data as a character
#' data.frame
#'
#' Row indices represent feature, columns represent feature attributions.
#'
#' \describe{
#'
#'    \item{.Data}{data-frame data, inherited from the data.frame class.}
#'
#'    \item{row.names}{
#'	     Also inherited from the data.frame class;
#'       it should contain the sample names.
#'    }
#'
#'    \item{names}{Inherited from the data.frame class.}
#'
#'  }
#'
#' @name Feature_table-class
#' @rdname Feature_table-class
#' @exportClass Feature_table
setClass("Feature_table", contains = "data.frame")
################################################################################
# Use setClassUnion to define the unholy NULL-data union as a virtual class.
# This is a way of dealing with the expected scenarios in which one or more of
# the component data classes is not available, in which case NULL will be used
# instead.
################################################################################
#' @keywords internal
setClassUnion("Profile_tableOrNULL", c("Profile_table", "NULL"))
#' @keywords internal
setClassUnion("Sample_dataOrNULL", c("Sample_data", "NULL"))
#' @keywords internal
setClassUnion("Feature_tableOrNULL", c("Feature_table", "NULL"))
################################################################################
#' The main experiment-level class for MyDataSet data
#'
#' Contains all currently-supported component data classes:
#' \code{\link{Profile_table-class}},
#' \code{\link{Sample_data-class}},
#' \code{\link{Feature_table-class}},
#' \code{MyDataSet-class} object the only data argument required for analysis and plotting
#' functions (although there are many options and parameter arguments available
#' to you).
#'
#' In the case of missing component data, the slots are set to \code{NULL}. As
#' soon as a \code{MyDataSet-class} object is to be updated with new component
#' data (previously missing/\code{NULL} or not), the indices of all components
#' are re-checked for compatibility and trimmed if necessary. This is to ensure
#' by design that components describe the same features/samples, and also that these
#' trimming/validity checks do not need to be repeated in downstream analyses.
#'
#' slots:
#' \describe{
#'    \item{Profile_table}{a single object of class Profile_table}
#'    \item{Sample_data}{ a single object of class Sample_data}
#'    \item{Feature_table}{ a single object of class Feature_table}
#' }
#'
#' @name MyDataSet-class
#' @rdname MyDataSet-class
#' @exportClass MyDataSet
#'
setClass(
    Class = "MyDataSet",
    representation = representation(
        Profile_table = "Profile_tableOrNULL",
        Sample_data = "Sample_dataOrNULL",
        Feature_table = "Feature_tableOrNULL"
    ),
    prototype = prototype(
        Profile_table = NULL,
        Sample_data = NULL,
        Feature_table = NULL
    )
)
################################################################################
#' Build MyDataSet-class objects
#'
#' This the constructor to build the [`MyDataSet-class`] object.
#' @param Profile_table, Numeric matrix; (Required)a Matrix of expression data, whose Row is FeatureID and Column is SampleID..
#' @param Sample_data, Data.frame; (Required)a dataframe. of Metadata(1st column must be "SampleID"), containing Group information and also environmental factors(biological factors).
#' @param Feature_table, Data.frame; the feature of Profile.
#'
#' @seealso [phyloseq::phyloseq()]
#' @name MyDataSet
#' @export
#' @return  a [`MyDataSet-class`] object.
#' @examples
#' \donttest{
#' library(dplyr)
#' Profile <- data.table::fread(system.file("extdata", "Species_relative_abundance.tsv", package="MyRtools"))  %>% tibble::column_to_rownames("V1")
#' Metadata <- read.csv(system.file("extdata", "Metadata.csv", package="MyRtools"))
#' Feature <- read.csv(system.file("extdata", "Species_feature.csv", package="MyRtools")) %>% tibble::column_to_rownames("Species")
#'
#' PROF <- Profile_table(as.matrix(Profile), feature_are_row=TRUE)
#' META <- Sample_data(Metadata)
#' FEATURE <- Feature_table(Feature)
#'
#' mydataset <- MyDataSet(Profile_table=PROF,
#'                        Sample_data=META,
#'                        Feature_table=FEATURE)
#' mydataset
#' }
#'
MyDataSet <- function(Profile_table = NULL,
                      Sample_data = NULL,
                      Feature_table = NULL){

    if(is.null(Profile_table)){
        stop("Profile_table is required")
    }

    if(is.null(Sample_data)){
        message("Sample_data is missing please check your input if you need this data")
    }

    if(is.null(Feature_table)){
        message("Feature_table is missing please check your input if you need this data")
    }

    new(
        "MyDataSet",
        Profile_table = Profile_table,
        Sample_data = Sample_data,
        Feature_table = Feature_table
    )
}
################################################################################
