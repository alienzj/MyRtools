################################################################################
# https://github.com/joey711/MyDataSet/blob/master/R/
################################################################################
#' Build or access the Profile_table
#'
#' This is the suggested method for both constructing and accessing
#' Operational Feature profile (\code{\link{Profile_table-class}}) objects.
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
    # Want dummy feature/sample index names if missing
    if(feature_are_rows){
        if(is.null(rownames(prftab))){
            rownames(prftab) <- paste("Feature", 1:nrow(prftab), sep="")
        }
        if(is.null(colnames(prftab))){
            colnames(prftab) <- paste("Sample", 1:ncol(prftab), sep="")
        }
    } else {
        if(is.null(rownames(prftab))){
            rownames(prftab) <- paste("Sample",1:nrow(prftab),sep="")
        }
        if(is.null(colnames(prftab))){
            colnames(prftab) <- paste("Feature",1:ncol(prftab),sep="")
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
#' Returns the total number of individuals observed from each feature.
#'
#' A convenience function equivalent to rowSums or colSums, but where
#' the orientation of the Profile_table is automatically handled.
#'
#' @usage feature_sums(x)
#'
#' @param x \code{\link{Profile_table-class}}, or \code{\link{MyDataSet-class}}.
#'
#' @return A \code{\link{numeric-class}} with length equal to the number of species
#'  in the table, name indicated the feature ID, and value equal to the sum of
#'  all individuals observed for each feature in \code{x}.
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
#' character matrix representing the description of
#' feature in the experiment.
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
# Constructor; for creating Feature_table from a matrix.
#' @rdname Feature_table-methods
#' @aliases Feature_table,matrix-method
setMethod("Feature_table", "matrix", function(object){
    # Want dummy species/taxa index names if missing
    if(is.null(rownames(object))){
        rownames(object) <- paste("Feature", 1:nrow(object), sep="")
    }
    if(is.null(colnames(object))){
        colnames(object) <- paste("Attribution", 1:ncol(object), sep="")
    }
    # instantiate as Feature_table
    return(new("Feature_table", object))
})
# Constructor; coerce to matrix, then pass on for creating Feature_table
#' @rdname Feature_table-methods
#' @aliases Feature_table,data.frame-method
setMethod("Feature_table", "data.frame", function(object){
    return(new("Feature_table", object))
})
################################################################################
#' Access feature_are_rows slot from Profile_table objects.
#'
#' @usage feature_are_rows(mydataset)
#'
#' @param mydataset (Required). \code{\link{MyDataSet-class}}, or \code{\link{Profile_table-class}}.
#'
#' @return A logical indicating the orientation of the Profile_table
#'
#' @seealso \code{\link{Profile_table}}
#' @rdname feature_are_rows-methods
#' @docType methods
#' @export
#' @aliases feature_are_rows feature_are_rows
setGeneric("feature_are_rows", function(mydataset) standardGeneric("feature_are_rows"))
#' @rdname feature_are_rows-methods
#' @aliases feature_are_rows,ANY-method
setMethod("feature_are_rows", "ANY", function(mydataset){NULL})
#' @rdname feature_are_rows-methods
#' @aliases feature_are_rows,Profile_table-method
setMethod("feature_are_rows", "Profile_table", function(mydataset){mydataset@feature_are_rows})
#' @rdname feature_are_rows-methods
#' @aliases feature_are_rows,MyDataSet-method
setMethod("feature_are_rows", "MyDataSet", function(mydataset){
    feature_are_rows(Profile_table(mydataset))
})
################################################################################
#' Get feature names.
#'
#' @usage feature_names(mydataset)
#'
#' @param mydataset \code{\link{MyDataSet-class}}, \code{\link{Profile_table-class}},
#'  \code{\link{Feature_table-class}}
#'
#' @return A character vector of the names of the species in \code{mydataset}.
#'
#'
#' @rdname feature_names-methods
#' @docType methods
#' @export
#'
#' @examples
#'
setGeneric("feature_names", function(mydataset) standardGeneric("feature_names"))
#' @rdname feature_names-methods
#' @aliases feature_names,ANY-method
setMethod("feature_names", "ANY", function(mydataset){ return(NULL) })
#' @rdname feature_names-methods
#' @aliases feature_names,MyDataSet-method
setMethod("feature_names", "MyDataSet", function(mydataset){
    feature_names(Profile_table(mydataset))
})
#' @rdname feature_names-methods
#' @aliases feature_names,Profile_table-method
setMethod("feature_names", "Profile_table", function(mydataset){
    if( feature_are_rows(mydataset) ){
        return( rownames(mydataset) )
    } else {
        return( colnames(mydataset) )
    }
})
#' @rdname feature_names-methods
#' @aliases feature_names,Feature_table-method
setMethod("feature_names", "Feature_table", function(mydataset) rownames(mydataset) )
#' @rdname feature_names-methods
#' @aliases feature_names,Sample_data-method
setMethod("feature_names", "Sample_data", function(mydataset) NULL )
################################################################################
#' Get the number of samples.
#'
#' @usage nsamples(mydataset)
#'
#' @param mydataset A \code{\link{MyDataSet-class}}, \code{\link{Sample_data}},
#'  or \code{\link{Profile_table-class}}.
#'
#' @return An integer indicating the total number of samples.
#'
#' @seealso \code{\link{feature_names}}, \code{\link{sample_names}},
#'  \code{\link{nfeature}}
#'
#' @rdname nsamples-methods
#' @docType methods
#' @export
#'
#' @examples
#'
setGeneric("nsamples", function(mydataset) standardGeneric("nsamples"))
#' @rdname nsamples-methods
#' @aliases nsamples,ANY-method
setMethod("nsamples", "ANY", function(mydataset){ return(NULL) })
#' @rdname nsamples-methods
#' @aliases nsamples,MyDataSet-method
setMethod("nsamples", "MyDataSet", function(mydataset){
    # dispatch to core, required component, Profile_table
    nsamples(Profile_table(mydataset))
})
#' @rdname nsamples-methods
#' @aliases nsamples,Profile_table-method
setMethod("nsamples", "Profile_table", function(mydataset){
    if( feature_are_rows(mydataset) ){
        return( ncol(mydataset) )
    } else {
        return( nrow(mydataset) )
    }
})
#' @rdname nsamples-methods
#' @aliases nsamples,Sample_data-method
setMethod("nsamples", "Sample_data", function(mydataset) nrow(mydataset) )
################################################################################
################################################################################
#' Get sample names.
#'
#' @usage sample_names(mydataset)
#'
#' @param mydataset (Required). A \code{\link{MyDataSet-class}}, \code{\link{Sample_data}},
#'  or \code{\link{Profile_table-class}}.
#'
#' @return A character vector. The names of the samples in \code{mydataset}.
#'
#' @aliases sample_names
#'
#' @rdname sample_names-methods
#' @docType methods
#' @export
#'
#' @examples
#'
setGeneric("sample_names", function(mydataset) standardGeneric("sample_names"))
# Unless otherwise specified, this should return a value of NULL
# That way, objects that do not explicitly describe samples all
# behave in the same (returning NULL) way.
#' @rdname sample_names-methods
#' @aliases sample_names,ANY-method
setMethod("sample_names", "ANY", function(mydataset){ return(NULL) })
#' @rdname sample_names-methods
#' @aliases sample_names,MyDataSet-method
setMethod("sample_names", "MyDataSet", function(mydataset){
    sample_names(Profile_table(mydataset))
})
#' @rdname sample_names-methods
#' @aliases sample_names,Sample_data-method
setMethod("sample_names", "Sample_data", function(mydataset) rownames(mydataset) )
#' @rdname sample_names-methods
#' @aliases sample_names,Profile_table-method
setMethod("sample_names", "Profile_table", function(mydataset){
    if( feature_are_rows(mydataset) ){
        return( colnames(mydataset) )
    } else {
        return( rownames(mydataset) )
    }
})
################################################################################
#' Get the number of feature.
#'
#' @usage nfeature(mydataset)
#'
#' @param mydataset \code{\link{MyDataSet-class}}, \code{\link{Profile_table-class}},
#'  \code{\link{Feature_table-class}}
#'
#' @return An integer indicating the number of feature.
#'
#' @seealso feature_names
#'
#' @rdname nfeature-methods
#' @docType methods
#' @export
#'
#' @examples
#'
setGeneric("nfeature", function(mydataset) standardGeneric("nfeature"))
#' @rdname nfeature-methods
#' @aliases nfeature,ANY-method
setMethod("nfeature", "ANY", function(mydataset){ return(NULL) })
#' @rdname nfeature-methods
#' @aliases nfeature,MyDataSet-method
setMethod("nfeature", "MyDataSet", function(mydataset){
    nfeature(Profile_table(mydataset))
})
#' @rdname nfeature-methods
#' @aliases nfeature,Profile_table-method
setMethod("nfeature", "Profile_table", function(mydataset){
    if( feature_are_rows(mydataset) ){
        return( nrow(mydataset) )
    } else {
        return( ncol(mydataset) )
    }
})
#' @rdname nfeature-methods
#' @aliases nfeature,Feature_table-method
setMethod("nfeature", "Feature_table", function(mydataset){
    nrow(mydataset)
})
################################################################################
#' @rdname show-methods
setMethod("show", "Profile_table", function(object){
    # print Profile_table (always there).
    cat(paste("Profile Table:          [", nfeature(object), " feature and ",
              nsamples(object), " samples]", sep = ""), fill = TRUE)
    if( feature_are_rows(object) ){
        cat("                     feature are rows", fill=TRUE)
    } else {
        cat("                     feature are columns", fill=TRUE)
    }
    show(as(object, "matrix"))
})
############################################################################
#' @rdname show-methods
setMethod("show", "Sample_data", function(object){
    cat(paste("Sample Data:        [", dim(Sample_data(object))[1], " samples by ",
              dim(Sample_data(object))[2],
              " sample variables]:", sep = ""),
        fill = TRUE)
    show(as(object, "data.frame"))
})
############################################################################
#' @rdname show-methods
setMethod("show", "Feature_Table", function(object){
    cat(paste("Feature Table:     [", dim(object)[1], " Feature by ",
              dim(object)[2],
              " Feature ranks]:", sep = ""),
        fill = TRUE)
    show(as(object, "data.frame"))
})
############################################################################
#' method extensions to show for MyDataSet objects.
#'
#' See the general documentation of \code{\link[methods]{show}} method for
#' expected behavior.
#'
#' @seealso \code{\link[methods]{show}}
#'
#' @inheritParams methods::show
#' @export
#' @rdname show-methods
#' @examples
#'
setMethod("show", "MyDataSet", function(object){
    cat("MyDataSet-class experiment-level object", fill=TRUE)
    # print Profile_table (always there).
    cat(paste("Profile_table()   Profile Table:         [ ", nfeature(Profile_table(object)), " Features and ",
              nsamples(Profile_table(object)), " Samples ]", sep = ""), fill = TRUE)

    # print Sample Data if there
    if(!is.null(Sample_data(object, FALSE))){
        cat(paste("Sample_data() Sample Data:       [ ", dim(Sample_data(object))[1], " Samples by ",
                  dim(Sample_data(object))[2],
                  " Sample variables ]", sep = ""), fill = TRUE)
    }

    # print Feature table if there
    if(!is.null(Feature_table(object, FALSE))){
        cat(paste("Feature_table()   Feature Table:    [ ", dim(Feature_table(object))[1], " Features by ",
                  dim(Feature_table(object))[2],
                  " Feature Attributions ]", sep = ""), fill = TRUE)
    }
})
############################################################################
