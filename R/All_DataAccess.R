################################################################################
# https://github.com/joey711/MyDataSet/blob/master/R/
################################################################################
#' Build or access the Profile_table
#'
#' This is the suggested method for both constructing and accessing
#' Operational Feature abundance (\code{\link{Profile_table-class}}) objects.
#' When the first
#' argument is a matrix, Profile_table() will attempt to create and return an
#' Profile_table-class object,
#' which further depends on whether or not \code{feature_are_rows} is provided as an
#' additional argument.
#' Alternatively, if the first argument is an experiment-level (\code{\link{MyDataSet-class}})
#' object, then the corresponding \code{Profile_table} is returned.
#'
#' @usage Profile_table(object, feature_are_rows, errorIfNULL=TRUE)
#'
#' @param object (Required). An integer matrix, \code{\link{Profile_table-class}},
#'  or \code{\link{MyDataSet-class}}.
#'
#' @param feature_are_rows (Conditionally optional). Logical; of length 1. Ignored
#'  unless \code{object} is a matrix, in which case it is is required.
#'
#' @param errorIfNULL (Optional). Logical. Should the accessor stop with
#'  an error if the slot is empty (\code{NULL})? Default \code{TRUE}. Ignored
#'  if \code{object} argument is a matrix (constructor invoked instead).
#'
#' @return An \code{\link{Profile_table-class}} object.
#'
#' @docType methods
#' @rdname Profile_table-methods
#' @export
#' @examples
#'
#'
setGeneric("Profile_table", function(object, feature_are_rows, errorIfNULL=TRUE){
    standardGeneric("Profile_table")
})
# Access the Profile_table slot.
#' @aliases Profile_table,MyDataSet-method
#' @rdname Profile_table-methods
setMethod("Profile_table", "MyDataSet", function(object, errorIfNULL=TRUE){
    access(object, "Profile_table", errorIfNULL)
})
# return the Profile_table as-is.
#' @aliases Profile_table,Profile_table-method
#' @rdname Profile_table-methods
setMethod("Profile_table", "Profile_table", function(object, errorIfNULL=TRUE){ return(object) })
# Instantiate an Profile_table from a raw abundance matrix.
#' @aliases Profile_table,matrix-method
#' @rdname Profile_table-methods
setMethod("Profile_table", "matrix", function(object, feature_are_rows){
    # instantiate first to check validity
    prftab <- new("Profile_table", object, feature_are_rows=feature_are_rows)
    # Want dummy species/sample index names if missing
    if(feature_are_rows){
        if(is.null(rownames(prftab))){
            rownames(prftab) <- paste("ft", 1:nrow(prftab), sep="")
        }
        if(is.null(colnames(prftab))){
            colnames(prftab) <- paste("sa", 1:ncol(prftab), sep="")
        }
    } else {
        if(is.null(rownames(prftab))){
            rownames(prftab) <- paste("sa",1:nrow(prftab),sep="")
        }
        if(is.null(colnames(prftab))){
            colnames(prftab) <- paste("ft",1:ncol(prftab),sep="")
        }
    }
    return(prftab)
})
# # # Convert to matrix, then dispatch.
#' @aliases Profile_table,data.frame-method
#' @rdname Profile_table-methods
setMethod("Profile_table", "data.frame", function(object, feature_are_rows){
    Profile_table(as(object, "matrix"), feature_are_rows)
})
# Any less-specific class, not inherited by those above.
#' @aliases Profile_table,ANY-method
#' @rdname Profile_table-methods
setMethod("Profile_table", "ANY", function(object, errorIfNULL=TRUE){
    access(object, "Profile_table", errorIfNULL)
})
################################################################################
#' Returns the total number of individuals observed from each species/taxa/OTU.
#'
#' A convenience function equivalent to rowSums or colSums, but where
#' the orientation of the Profile_table is automatically handled.
#'
#' @usage feature_sums(x)
#'
#' @param x \code{\link{Profile_table-class}}, or \code{\link{MyDataSet-class}}.
#'
#' @return A \code{\link{numeric-class}} with length equal to the number of species
#'  in the table, name indicated the taxa ID, and value equal to the sum of
#'  all individuals observed for each taxa in \code{x}.
#'
#' @export
#' @examples
#'
feature_sums <- function(x){
    x <- Profile_table(x)
    if( feature_are_rows(x) ){
        rowSums(x)
    } else {
        colSums(x)
    }
}
################################################################################
#' Returns the total number of individuals observed from each sample.
#'
#' A convenience function equivalent to rowSums or colSums, but where
#' the orientation of the Profile_table is automatically handled.
#'
#' @usage sample_sums(x)
#'
#' @param x \code{\link{Profile_table-class}}, or \code{\link{MyDataSet-class}}.
#'
#' @return A named \code{\link{numeric-class}}
#'  length equal to the number of samples
#'  in the \code{x}, name indicating the sample ID, and value equal to the sum of
#'  all individuals observed for each sample in \code{x}.
#'
#' @export
#' @examples
#'
sample_sums <- function(x){
    x <- Profile_table(x)
    if( feature_are_rows(x) ){
        colSums(x)
    } else {
        rowSums(x)
    }
}
################################################################################
#' Build or access Sample_data.
#'
#' This is the suggested method for both constructing and accessing a table
#' of sample-level variables (\code{\link{Sample_data-class}}),
#' which in the \code{\link{MyDataSet-package}} is represented as a special
#' extension of the \code{\link{data.frame-class}}.
#' When the
#' argument is a \code{\link{data.frame}}, \code{Sample_data} will create
#' a Sample_data-class object.
#' In this case, the rows should be named to match the
#' \code{\link{sample_names}} of the other objects to which it will ultimately be paired.
#' Alternatively, if the first argument is an experiment-level (\code{\link{MyDataSet-class}})
#' object, then the corresponding \code{Sample_data} is returned.
#' Like other accessors (see See Also, below), the default behavior of this method
#' is to stop with an
#' error if \code{object} is a \code{MyDataSet-class} but does not
#' contain a \code{Sample_data}.
#'
#' @usage Sample_data(object, errorIfNULL=TRUE)
#'
#' @param object (Required). A \code{\link{data.frame-class}},
#'  or a \code{\link{MyDataSet-class}} object.
#'
#' @param errorIfNULL (Optional). Logical. Should the accessor stop with
#'  an error if the slot is empty (\code{NULL})? Default \code{TRUE}.
#'
#' @return A \code{\link{Sample_data-class}} object
#' representing the sample variates of an experiment.
#'
#' @aliases Sample_data
#'
#' @rdname Sample_data-methods
#' @docType methods
#' @export
#'
#' @examples
#'
setGeneric("Sample_data", function(object, errorIfNULL=TRUE) standardGeneric("Sample_data"))
#' @rdname Sample_data-methods
#' @aliases Sample_data,ANY-method
setMethod("Sample_data", "ANY", function(object, errorIfNULL=TRUE){
    access(object, "Sample_data", errorIfNULL)
})
# constructor; for creating Sample_data from a data.frame
#' @rdname Sample_data-methods
#' @aliases Sample_data,data.frame-method
setMethod("Sample_data", "data.frame", function(object){
    # Make sure there are no phantom levels in categorical variables
    object <- reconcile_categories(object)

    # instantiate first to check validity
    SM <- new("Sample_data", object)
    return(SM)
})
################################################################################
#' Cleans absent levels in Sample_data/data.frame.
#'
#' This is used internally by the builder method, \code{\link{Sample_data}}, to
#' ensure that the factors describing categorical variables in a data.frame or
#' Sample_data object are free of extra levels that can plague downstream plots
#' analysis.
#'
#' @usage reconcile_categories(DFSM)
#'
#' @param DFSM (Required). A \code{data.frame} or \code{Sample_data} object that needs to be cleaned.
#'
#' @return A single \code{data.frame} object. Even if the input argument is a \code{Sample_data},
#'  the return is a \code{data.frame}. Because this is intended to be used internally by
#'  the builder method, it cannot also call the builder function to re-build
#'  the cleaned \code{Sample_data}.
#'
#' @keywords internal
#'
#' @examples
#'
reconcile_categories <- function(DFSM){
    DF = as(DFSM, "data.frame")
    return(droplevels(DF))
}
################################################################################
#' Build or access the Feature_table.
#'
#' This is the suggested method for both constructing and accessing a table of
#' Feature names, organized with ranks as columns (\code{\link{Feature_table-class}}).
#' When the argument is a character matrix, Feature_table() will create and return a
#' \code{\link{Feature_table-class}} object.
#' Alternatively, if the first argument is an experiment-level (\code{\link{MyDataSet-class}})
#' object, then the corresponding \code{Feature_table} is returned.
#' Like other accessors (see See Also, below), the default behavior of this method
#' is to stop with an
#' error if \code{object} is a \code{MyDataSet-class} but does not
#' contain a \code{Feature_table}.
#'
#' @usage Feature_table(object, errorIfNULL=TRUE)
#'
#' @param object An object among the set of classes defined by the MyDataSet
#' package that contain Feature_table
#'
#' @param errorIfNULL (Optional). Logical. Should the accessor stop with
#'  an error if the slot is empty (\code{NULL})? Default \code{TRUE}.
#'
#' @return A \code{\link{Feature_table-class}} object.
#' It is either grabbed from the relevant slot
#' if \code{object} is complex, or built anew if \code{object} is a
#' character matrix representing the taxonomic classification of
#' species in the experiment.
#'
#' @rdname Feature_table-methods
#' @docType methods
#' @export
#'
#' @examples
#'
setGeneric("Feature_table", function(object, errorIfNULL=TRUE) standardGeneric("Feature_table"))
#' @rdname Feature_table-methods
#' @aliases Feature_table,ANY-method
setMethod("Feature_table",  "ANY", function(object, errorIfNULL=TRUE){
    access(object, "Feature_table", errorIfNULL)
})
# Constructor; coerce to matrix, then pass on for creating Feature_table
#' @rdname Feature_table-methods
#' @aliases Feature_table,data.frame-method
setMethod("Feature_table", "data.frame", function(object){
    return(new("Feature_table", object))
})
################################################################################
#' Show the component objects classes and slot names.
#'
#' There are no arguments to this function. It returns a named character
#' when called, which can then be used for tests of component data types, etc.
#'
#' @usage get.component.classes()
#'
#' @return a character vector of the component objects classes, where each
#' element is named by the corresponding slot name in the MyDataSet-class.
#'
#' @keywords internal
#'
#' @examples #
#'
get.component.classes <- function(){
    # define classes vector
    component.classes <- c("Profile_table", "Sample_data", "Feature_table")
    # the names of component.classes needs to be the slot names to match getSlots / splat
    names(component.classes) <- c("Profile_table", "Sample_data", "Feature_table")
    return(component.classes)
}
# Explicitly define components/slots that describe taxa.
#' @keywords internal
feature.components <- function(){
    # define classes vector
    component.classes <- c("Profile_table", "Feature_table")
    # the names of component.classes needs to be the slot names to match getSlots / splat
    names(component.classes) <- c("Profile_table", "Feature_table")
    return(component.classes)
}
# Explicitly define components/slots that describe samples.
#' @keywords internal
sample.components <- function(){
    # define classes vector
    component.classes <- c("Profile_table", "Sample_data")
    # the names of component.classes needs to be the slot names to match getSlots / splat
    names(component.classes) <- c("Profile_table", "Sample_data")
    return(component.classes)
}
# Returns TRUE if x is a component class, FALSE otherwise.
# This shows up over and over again in data infrastructure
#' @keywords internal
is.component.class <- function(x){
    inherits(x, get.component.classes())
}
################################################################################
#' Universal slot accessor function for MyDataSet-class.
#'
#' This function is used internally by many accessors and in
#' many functions/methods that need to access a particular type of component data.
#' If something is wrong, or the slot is missing, the expected behavior is that
#' this function will return NULL. Thus, the output can be tested by
#' \code{\link{is.null}} as verification of the presence of a particular
#' data component.
#' the default behavior is not to stop with an error if the desired slot is empty.
#' In all cases this is controlled by the \code{errorIfNULL} argument, which can
#' be set to \code{TRUE} if an error is desired.
#'
#' @usage access(MyDataSet, slot, errorIfNULL=FALSE)
#'
#' @param MyDataSet (Required). \code{\link{MyDataSet-class}}.
#'
#' @param slot (Required). A character string indicating the slot (not data class)
#'  of the component data type that is desired.
#'
#' @param errorIfNULL (Optional). Logical. Should the accessor stop with
#'  an error if the slot is empty (\code{NULL})? Default \code{FALSE}.
#'
#' @return Returns the component object specified by the argument \code{slot}.
#'  Returns NULL if slot does not exist. Returns \code{MyDataSet} as-is
#'  if it is a component class that already matches the slot name.
#'
#' @export
#' @examples
#'
access <- function(mydataset, slot, errorIfNULL=FALSE){
    if( is.component.class(mydataset) ){
        # If mydataset is a component class, might return as-is. Depends on slot.
        if( inherits(mydataset, get.component.classes()[slot]) ){
            # if slot-name matches, return mydataset as-is.
            out = mydataset
        } else {
            # If slot/component mismatch, set out to NULL. Test later if this is an error.
            out = NULL
        }
    } else if(!slot %in% slotNames(mydataset) ){
        # If slot is invalid, set out to NULL. Test later if this is an error.
        out = NULL
    } else {
        # By elimination, must be valid. Access slot
        out = eval(parse(text=paste("mydataset@", slot, sep="")))
    }
    if( errorIfNULL & is.null(out) ){
        # Only error regarding a NULL return value if errorIfNULL is TRUE.
        stop(slot, " slot is empty.")
    }
    return(out)
}
################################################################################
