################################################################################
# https://github.com/joey711/MyDataSet/blob/master/R/
################################################################################
#' Assign a new Profile Table to \code{x}
#'
#' @usage Profile_table(x) <- value
#'
#' @param x (Required). \code{\link{MyDataSet-class}}
#' @param value (Required).
#'  \code{\link{Profile_table-class}}
#'  or
#'  \code{\link{MyDataSet-class}}.
#'
#' @export
#' @docType methods
#' @rdname assign-Profile_table
#' @aliases assign-Profile_table
#'
#' @examples
#'
setGeneric("Profile_table<-", function(x, value) standardGeneric("Profile_table<-"))
#' @rdname assign-Profile_table
#' @aliases Profile_table<-,MyDataSet,Profile_table-method
setMethod("Profile_table<-", c("MyDataSet", "Profile_table"), function(x, value){
    MyDataSet(value, x@Sample_data, x@Feature_table)
})
#' @rdname assign-Profile_table
#' @aliases Profile_table<-,Profile_table,Profile_table-method
setMethod("Profile_table<-", c("Profile_table", "Profile_table"), function(x, value){ value })
#' @rdname assign-Profile_table
#' @aliases Profile_table<-,MyDataSet,MyDataSet-method
setMethod("Profile_table<-", c("MyDataSet", "MyDataSet"), function(x, value){
    MyDataSet(Profile_table(value), x@Sample_data, x@Feature_table)
})
################################################################################
#' Manually change feature_are_rows through assignment.
#'
#' The feature_are_rows slot is a logical indicating the orientation of the
#' abundance table contained in object \code{x}.
#'
#' @usage feature_are_rows(x) <- value
#'
#' @param x \code{\link{Profile_table-class}} or \code{\link{MyDataSet-class}}
#'
#' @param value A logical of length equal to 1. If \code{length(value) > 1},
#'  the additional elements will be ignored. Only the first element is assigned
#'  to the feature_are_rows slot.
#'
#' @export
#' @docType methods
#' @rdname assign-feature_are_rows
#' @aliases assign-feature_are_rows feature_are_rows<-
#'
#' @examples
#'
setGeneric("feature_are_rows<-", function(x, value){
    standardGeneric("feature_are_rows<-")
})
#' @rdname assign-feature_are_rows
#' @aliases feature_are_rows<-,Profile_table,logical-method
setMethod("feature_are_rows<-", c("Profile_table", "logical"), function(x, value){
    x@feature_are_rows <- value[1]
    return(x)
})
#' @rdname assign-feature_are_rows
#' @aliases feature_are_rows<-,MyDataSet,logical-method
setMethod("feature_are_rows<-", c("MyDataSet", "logical"), function(x, value){
    feature_are_rows(Profile_table(x)) <- value
    return(x)
})
################################################################################
#' Assign (new) Sample_data to \code{x}
#'
#' This replaces the current \code{Sample_data} component of \code{x} with
#' \code{value}, if \code{value} is a \code{\link{Sample_data-class}}. However,
#' if \code{value} is a \code{data.frame}, then \code{value} is first coerced to
#' a \code{\link{Sample_data-class}}, and then assigned. Alternatively, if
#' \code{value} is \code{\link{MyDataSet-class}}, then the
#' \code{\link{Sample_data}} component will first be accessed from \code{value}
#'  and then assigned. This makes possible some concise assignment/replacement
#'  statements when adjusting, modifying, or building subsets of
#'  experiment-level data. See some examples below.
#'
#' Internally, this re-builds the \code{\link{MyDataSet-class}} object using
#' the standard \code{\link{MyDataSet}} constructor. Thus, index mismatches
#' between sample-describing components will not be allowed, and subsetting
#' will occurr automatically such that only the intersection of sample IDs
#' are included in any components. This has the added benefit of re-checking
#' (internally) for any other issues.
#'
#' @usage Sample_data(x) <- value
#'
#' @param x (Required). \code{\link{MyDataSet-class}}. The object to modify.
#' @param value (Required). Either a \code{\link{Sample_data-class}},
#'  a \code{data.frame} that can be coerced into \code{\link{Sample_data-class}},
#'  or a \code{\link{MyDataSet-class}} that contains a
#'  suitable \code{Sample_data} component to assign to \code{x}. If unsure,
#'  try \code{\link{Sample_data}}\code{(value)}, which should return a
#'  \code{\link{Sample_data-class}} object without error.
#'
#' @return No return. This is an assignment statement.
#'
#' @export
#' @rdname assign-Sample_data
#' @aliases assign-Sample_data Sample_data<-
#' @examples
#'
"Sample_data<-" <- function(x, value){
    if( !inherits(value, "Sample_data") ){
        value <- Sample_data(value)
    }
    MyDataSet(x@Profile_table, value, x@Feature_table)
}
################################################################################
#' Assign a (new) Feature Table to \code{x}
#'
#' @usage Feature_table(x) <- value
#'
#' @param x (Required). \code{\link{MyDataSet-class}}
#' @param value (Required). \code{\link{Feature_table-class}}.
#'
#' @export
#' @rdname assign-Feature_table
#' @aliases assign-Feature_table Feature_table<-
#' @examples
#'
setGeneric("Feature_table<-", function(x, value) standardGeneric("Feature_table<-"))
#' @rdname assign-Feature_table
#' @aliases Feature_table<-,MyDataSet,Feature_table-method
setMethod("Feature_table<-", c("MyDataSet", "Feature_table"), function(x, value){
    MyDataSet(x@Profile_table, x@Sample_data, value)
})
#' @rdname assign-Feature_table
#' @aliases Feature_table<-,MyDataSet,ANY-method
setMethod("Feature_table<-", c("MyDataSet", "ANY"), function(x, value){
    MyDataSet(x@Profile_table, x@Sample_data, Feature_table(value, FALSE))
})
#' @rdname assign-Feature_table
#' @aliases Feature_table<-,Feature_table,Feature_table-method
setMethod("Feature_table<-", c("Feature_table", "Feature_table"), function(x, value){
    # Asign as-is.
    value
})
#' @rdname assign-Feature_table
#' @aliases Feature_table<-,Feature_table,ANY-method
setMethod("Feature_table<-", c("Feature_table", "ANY"), function(x, value){
    Feature_table(value, FALSE)
})
################################################################################
#' Replace Feature identifier names
#'
#' @usage feature_names(x) <- value
#'
#' @param x (Required). An object defined by the \code{\link{MyDataSet-package}}
#' 	that describes Features in some way.
#' @param value (Required). A character vector
#'  to replace the current \code{\link{feature_names}}.
#'
#' @export
#' @docType methods
#' @rdname assign-feature_names
#' @aliases assign-feature_names feature_names<-
#'
#' @examples
#'
setGeneric("feature_names<-", function(x, value){
    if( anyDuplicated(value) ){
        stop("feature_names<-: You are attempting to assign duplicated feature_names")
    }
    standardGeneric("feature_names<-")
})
# Attempt to coerce value to a character vector. Remaining methods will require it.
#' @rdname assign-feature_names
#' @aliases feature_names<-,ANY,ANY-method
setMethod("feature_names<-", c("ANY", "ANY"), function(x, value){
    feature_names(x) <- as(value, "character")
    return(x)
})
# value is now character, but no specific method for first argumet
# return x unchanged.
#' @rdname assign-feature_names
#' @aliases feature_names<-,ANY,character-method
setMethod("feature_names<-", c("ANY", "character"), function(x, value){
    return(x)
})
#' @rdname assign-feature_names
#' @aliases feature_names<-,Profile_table,character-method
setMethod("feature_names<-", c("Profile_table", "character"), function(x, value){
    if( taxa_are_rows(x) ){
        rownames(x) <- value
    } else {
        colnames(x) <- value
    }
    return(x)
})
#' @rdname assign-feature_names
#' @aliases feature_names<-,Feature_table,character-method
setMethod("feature_names<-", c("Feature_table", "character"), function(x, value){
    rownames(x) <- value
    return(x)
})
#' @rdname assign-feature_names
#' @aliases feature_names<-,MyDataSet,character-method
setMethod("feature_names<-", c("MyDataSet", "character"), function(x, value){
    # dispatch on components
    feature_names(x@Profile_table) <- value
    return(x)
})
################################################################################
#' Replace Sample identifier names
#'
#' @usage sample_names(x) <- value
#'
#' @param x (Required). An object defined by the \code{\link{MyDataSet-package}}
#' 	that describes Samples in some way.
#' @param value (Required). A character vector
#'  to replace the current \code{\link{sample_names}}.
#'
#' @export
#' @docType methods
#' @rdname assign-sample_names
#' @aliases assign-sample_names sample_names<-
#'
#' @examples
#'
setGeneric("sample_names<-", function(x, value){
    if( anyDuplicated(value) ){
        stop("sample_names<-: You are attempting to assign duplicated sample_names")
    }
    standardGeneric("sample_names<-")
})
# Attempt to coerce value to a character vector. Remaining methods will require it.
#' @rdname assign-sample_names
#' @aliases sample_names<-,ANY,ANY-method
setMethod("sample_names<-", c("ANY", "ANY"), function(x, value){
    sample_names(x) <- as(value, "character")
    return(x)
})
# value is now character, but no specific method for first argumet
# return x unchanged.
#' @rdname assign-sample_names
#' @aliases sample_names<-,ANY,character-method
setMethod("sample_names<-", c("ANY", "character"), function(x, value){
    return(x)
})
#' @rdname assign-sample_names
#' @aliases sample_names<-,Profile_table,character-method
setMethod("sample_names<-", c("Profile_table", "character"), function(x, value){
    if( taxa_are_rows(x) ){
        colnames(x) <- value
    } else {
        rownames(x) <- value
    }
    return(x)
})
#' @rdname assign-sample_names
#' @aliases sample_names<-,sample_data,character-method
setMethod("sample_names<-", c("Sample_data", "character"), function(x, value){
    rownames(x) <- value
    return(x)
})
#' @rdname assign-sample_names
#' @aliases sample_names<-,MyDataSet,character-method
setMethod("sample_names<-", c("MyDataSet", "character"), function(x, value){
    # dispatch on components
    sample_names(x@Profile_table) <- value
    sample_names(x@Sample_data)  <- value
    return(x)
})
################################################################################
