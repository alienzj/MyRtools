################################################################################
# https://github.com/joey711/MyDataSet/blob/master/R/allClasses.R
################################################################################
#' Build MyDataSet-class objects
#'
#' This the constructor to build the [`MyDataSet-class`] object.
#' @param Profile_table, Numeric matrix; (Required)a Matrix of expression data, whose Row is FeatureID and Column is SampleID..
#' @param Sample_data, Data.frame; (Required)a dataframe. of Metadata(1st column must be "SampleID"), containing Group information and also environmental factors(biological factors).
#' @param Feature_table, Data.frame; the feature of Profile.
#'
#' @seealso [MyDataSet::MyDataSet()]
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
#' FEATURE <- Feature_table(as.matrix(Feature))
#'
#' mydataset <- MyDataSet(Profile_table=PROF,
#'                        Sample_data=META,
#'                        Feature_table=FEATURE)
#' mydataset
#' }
#'
MyDataSet <- function(...){

    arglist <- list(...)

    # Remove names from arglist. Will replace them based on their class
    names(arglist) <- NULL

    # ignore all but component data classes.
    arglist  <- arglist[sapply(arglist, is.component.class)]

    # Make the name-replaced, splatted list
    splatlist <- sapply(arglist, splat.MyDataSet.objects)

    ####################
    ## Need to determine whether to
    # (A) instantiate a new raw/uncleaned MyDataSet object, or
    # (B) return a single component, or
    # (C) to stop with an error because of incorrect argument types.
    if( length(splatlist) > length(get.component.classes()) ){
        stop("Too many components provided\n")
    } else if( length(names(splatlist)) > length(unique(names(splatlist))) ){
        stop("Only one of each component type allowed.\n")
    } else if( length(splatlist) == 1){
        return(arglist[[1]])
    } else {
        # Instantiate the MyDataSet-class object, ps.
        ps <- do.call("new", c(list(Class="MyDataSet"), splatlist) )
    }

    ####################
    ## Reconcile the feature and sample index names between components
    ## in the newly-minted MyDataSet object
    shared_feature <- intersect_feature(ps)
    shared_samples <- intersect_samples(ps)

    if( length(shared_feature) < 1 ){
        stop("Problem with feature indices among those you provided.\n",
             "Check using intersect() and feature_names()\n"
        )
    }
    if( length(shared_samples) < 1 ){
        stop("Problem with sample indices among those you provided.\n",
             "Check using intersect() and feature_names()\n"
        )
    }

    # Start with Feature indices
    ps <- prune_feature(shared_feature, ps)

    # Verify there is more than one component
    # that describes samples before attempting to reconcile.
    ps <- prune_samples(shared_samples, ps)

    # Force Both samples and feature indices to be in the same order.
    ps <- index_reorder(ps, "Both")

    return(ps)
}
################################################################################
# A relatively fast way to access from MyDataSet object components
# f - function name as character string
# mydataset - a MyDataSet object (MyDataSet-class instance)
#' @keywords internal
f_comp_ps <- function(f, mydataset){
    sapply(names(getSlots("MyDataSet")), function(i, ps){
        eval(parse(text=paste(f, "(ps@", i, ")", sep="")))
    }, mydataset)
}
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
# Explicitly define components/slots that describe feature.
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
#' Convert \code{\link{MyDataSet-class}} into a named list of its non-empty components.
#'
#' This is used in internal handling functions, and one of its key features
#' is that the names in the returned-list match the slot-names, which is useful
#' for constructing calls with language-computing functions like \code{\link{do.call}}.
#' Another useful aspect is that it only returns the contents of non-empty slots.
#' In general, this should only be used by MyDataSet-package developers. Standard
#' users should not need or use this function, and should use the accessors and
#' other tools that leave the multi-component object in one piece.
#'
#' @usage splat.MyDataSet.objects(x)
#'
#' @param x A \code{\link{MyDataSet-class}} object. Alternatively, a component
#'  data object will work, resulting in named list of length 1.
#'
#' @return A named list, where each element is a component object that was contained
#' in the argument, \code{x}. Each element is named according to its slot-name in
#' the MyDataSet-object from which it is derived.
#' If \code{x} is already a component data object,
#' then a list of length (1) is returned, also named.
#'
#' @seealso merge_MyDataSet
#' @keywords internal
#' @examples #
splat.MyDataSet.objects <- function(x){
    if( is.component.class(x) ){
        # Check if class of x is among the component classes already (not MyDataSet-class)
        splatx <- list(x)
        names(splatx) <- names(which(sapply(get.component.classes(), function(cclass, x) inherits(x, cclass), x)))
    } else if( inherits(x, "MyDataSet") ){
        # Else, check if it inherits from MyDataSet, and if-so splat
        slotnames <- names(getSlots("MyDataSet"))
        allslots  <- sapply(slotnames, function(i, x){access(x, i, FALSE)}, x)
        splatx    <- allslots[!sapply(allslots, is.null)]
    } else {
        # Otherwise, who knows what it is, silently return NULL.
        return(NULL)
    }
    return(splatx)
}
################################################################################
#' Return the non-empty slot names of a MyDataSet object.
#'
#' Like \code{\link{getSlots}}, but returns the class name if argument
#' is component data object.
#'
#' @usage getslots.MyDataSet(mydataset)
#'
#' @param mydataset A \code{\link{MyDataSet-class}} object. If \code{mydataset} is a component
#'  data class, then just returns the class of \code{mydataset}.
#'
#' @return identical to getSlots. A named character vector of the slot classes
#' of a particular S4 class, where each element is named by the slot name it
#' represents. If \code{mydataset} is a component data object,
#' then a vector of length (1) is returned, named according to its slot name in
#' the \code{\link{MyDataSet-class}}.
#'
#' @seealso merge_MyDataSet
#' @export
#' @examples
#'
getslots.MyDataSet <- function(mydataset){
    names(splat.MyDataSet.objects(mydataset))
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
            out <- mydataset
        } else {
            # If slot/component mismatch, set out to NULL. Test later if this is an error.
            out <- NULL
        }
    } else if(!slot %in% slotNames(mydataset) ){
        # If slot is invalid, set out to NULL. Test later if this is an error.
        out <- NULL
    } else {
        # By elimination, must be valid. Access slot
        out <- eval(parse(text=paste("mydataset@", slot, sep="")))
    }
    if( errorIfNULL & is.null(out) ){
        # Only error regarding a NULL return value if errorIfNULL is TRUE.
        stop(slot, " slot is empty.")
    }
    return(out)
}
################################################################################
#' Returns the intersection of species and samples for the components of x
#'
#' This function is used internally as part of the infrastructure to ensure that
#' component data types in a MyDataSet-object have exactly the same feature.
#' It relies heavily on the \code{\link{Reduce}} function to determine the
#' strictly common species.
#'
#' @usage intersect_feature(x)
#'
#' @param x (Required). A \code{\link{MyDataSet-class}} object
#'  that contains 2 or more components
#'  that in-turn describe feature.
#'
#' @return Returns a character vector of only those feature that are present in
#'  all feature-describing components of \code{x}.
#'
#' @seealso \code{\link{Reduce}}, \code{\link{intersect}}
#' @keywords internal
#' @examples
#'
intersect_feature <- function(x){
    feature_vectors <- f_comp_ps("feature_names", x)
    feature_vectors <- feature_vectors[!sapply(feature_vectors, is.null)]
    return( Reduce("intersect", feature_vectors) )
}
#' @keywords internal
intersect_samples <- function(x){
    sample_vectors <- f_comp_ps("sample_names", x)
    sample_vectors <- sample_vectors[!sapply(sample_vectors, is.null)]
    return( Reduce("intersect", sample_vectors) )
}
################################################################################
#' Force index order of MyDataSet objects
#'
#' @usage index_reorder(ps, index_type)
#'
#' @param ps (Required). A \code{\link{MyDataSet-class}} instance.
#' @param index_type (Optional). A character string
#'  specifying the indices to properly order.
#'  Supported values are \code{c("Both", "Features", "Samples")}.
#'  Default is \code{"Both"}, meaning samples and Features indices
#'  will be checked/re-ordered.
#'
#' @keywords internal
#' @docType methods
#'
#' @examples
#' ## data("GlobalPatterns")
#' ## GP = index_reorder(GlobalPatterns)
setGeneric("index_reorder", function(ps, index_type) standardGeneric("index_reorder") )
#' @rdname index_reorder
#' @aliases index_reorder,MyDataSet-method
setMethod("index_reorder", "MyDataSet", function(ps, index_type="Both"){
    if( index_type %in% c("Both", "Feature") ){
        ## ENFORCE CONSISTENT ORDER OF FEATURE INDICES.
        torder <- feature_names(Profile_table(ps))
        if( !is.null(Feature_table(ps, FALSE)) ){
            ps@Feature_table <- Feature_table(Feature_table(ps)[torder, ])
        }
    }
    if( index_type %in% c("Both", "Samples") ){
        ## ENFORCE CONSISTENT ORDER OF SAMPLE INDICES
        if( !is.null(Sample_data(ps, FALSE)) ){
            # check first that ps has Sample_data
            if( !all(sample_names(Profile_table(ps)) == rownames(Sample_data(ps))) ){
                # Reorder the Sample_data rows so that they match the otu_table order.
                ps@Sample_data <- Sample_data(Sample_data(ps)[sample_names(Profile_table(ps)), ])
            }
        }
    }

    return(ps)
})
################################################################################
#' Define a subset of samples to keep in a MyDataSet object.
#'
#' An S4 Generic method for pruning/filtering unwanted samples
#' by defining those you want to keep.
#'
#' @usage prune_samples(samples, x)
#'
#' @param samples (Required). A character vector of the samples in object x that you want to
#' keep -- OR alternatively -- a logical vector where the kept samples are TRUE, and length
#' is equal to the number of samples in object x. If \code{samples} is a named
#' logical, the samples retained is based on those names. Make sure they are
#' compatible with the \code{sample_names} of the object you are modifying (\code{x}).
#'
#' @param x A MyDataSet object.
#'
#' @return The class of the object returned by \code{prune_samples} matches
#' the class of the MyDataSet object, \code{x}.
#'
#'
#' @rdname prune_samples-methods
#' @docType methods
#' @export
#' @examples
#'
setGeneric("prune_samples", function(samples, x) standardGeneric("prune_samples"))
#' @aliases prune_samples,character,Profile_table-method
#' @rdname prune_samples-methods
setMethod("prune_samples", signature("character", "Profile_table"), function(samples, x){
    if( setequal(samples, sample_names(x)) ){
        # If the sets of `samples` and sample_names are the same, return as-is.
        return(x)
    } else {
        samples <- intersect(samples, sample_names(x))
        if( feature_are_rows(x) ){
            return( x[, samples] )
        } else {
            return( x[samples, ] )
        }
    }
})
#' @aliases prune_samples,character,Sample_data-method
#' @rdname prune_samples-methods
setMethod("prune_samples", signature("character", "Sample_data"), function(samples, x){
    if( setequal(samples, sample_names(x)) ){
        # If the sets of `samples` and sample_names are the same, return as-is.
        return(x)
    } else {
        samples <- intersect(samples, sample_names(x))
        return(x[samples, ])
    }
})
#' @aliases prune_samples,character,MyDataSet-method
#' @rdname prune_samples-methods
setMethod("prune_samples", signature("character", "MyDataSet"), function(samples, x){
    # Re-define `samples` as the intersection of samples names for each component AND `samples`
    samples <- intersect(intersect_samples(x), samples)
    # Now prune each component.
    # All MyDataSet objects have an Profile_table slot, no need to test for existence.
    x@Profile_table <- Profile_table(prune_samples(samples, Profile_table(x)), feature_are_rows = TRUE)
    if( !is.null(x@Sample_data) ){
        # protect missing Sample_data component. Don't need to prune if empty
        x@Sample_data <- Sample_data(prune_samples(samples, Sample_data(x)))
    }
    # Force sample index order after pruning to be the same,
    # according to the same rules as in the constructor, MyDataSet()
    x <- index_reorder(x, index_type="Samples")
    return(x)
})
# A logical should specify the samples to keep, or not. Have same length as nsamples(x)
#' @aliases prune_samples,logical,ANY-method
#' @rdname prune_samples-methods
setMethod("prune_samples", signature("logical", "ANY"), function(samples, x){
    # Check that logical has same length as nsamples, stop if not.
    if( !identical(length(samples), nsamples(x)) ){
        stop("logical argument to samples is wrong length. Should equal nsamples(x)")
    } else {
        # Pass on to names-based prune_samples method
        return( prune_samples(sample_names(x)[samples], x) )
    }
})
################################################################################
#' Prune unwanted Features from a MyDataSet object.
#'
#' An S4 Generic method for removing (pruning) unwanted Features from MyDataSet
#' objectsas well as native MyDataSet package objects.
#' This is particularly useful for pruning a MyDataSet object that has
#' more than one component that describes Features.
#' Credit: the \code{phylo}-class version is adapted from
#' \href{http://cran.at.r-project.org/web/packages/picante/index.html}{prune.sample}.
#'
#' @usage prune_feature(feature, x)
#'
#' @param feature (Required). A character vector of the feature in object x that you want to
#' keep -- OR alternatively -- a logical vector where the kept feature are TRUE, and length
#' is equal to the number of feature in object x. If \code{feature} is a named
#' logical, the feature retained are based on those names. Make sure they are
#' compatible with the \code{feature_names} of the object you are modifying (\code{x}).
#'
#' @param x (Required). A MyDataSet object, including MyDataSet classes that represent feature.
#' If the function \code{\link{feature_names}} returns a non-\code{NULL} value, then your object
#' can be pruned by this function.
#'
#' @return The class of the object returned by \code{prune_feature} matches
#' the class of the argument, \code{x}.
#'
#' @seealso
#'
#'  \code{\link{prune_samples}}
#'
#'  \href{http://cran.at.r-project.org/web/packages/picante/index.html}{prune.sample}
#'
#' @rdname prune_feature-methods
#' @export
#' @examples
#'
setGeneric("prune_feature", function(feature, x) standardGeneric("prune_feature"))
#' @aliases prune_feature,NULL,ANY-method
#' @rdname prune_feature-methods
setMethod("prune_feature", signature("NULL", "ANY"), function(feature, x){
    return(x)
})
# Any prune_feature call w/ signature starting with a logical
# converts the logical to a character vector, and then dispatches
# to more specific method.
#' @aliases prune_feature,logical,ANY-method
#' @rdname prune_feature-methods
setMethod("prune_feature", signature("logical", "ANY"), function(feature, x){
    # Check that logical has same length as nfeature, stop if not.
    if( !identical(length(feature), nfeature(x)) ){
        stop("logical argument to feature is wrong length. Should equal nfeature(x)")
    } else {
        # Pass on to names-based prune_feature method
        return( prune_feature(feature_names(x)[feature], x) )
    }
})
#' @aliases prune_feature,character,Profile_table-method
#' @rdname prune_feature-methods
setMethod("prune_feature", signature("character", "Profile_table"), function(feature, x){
    if( setequal(feature, feature_names(x)) ){
        return(x)
    } else {
        feature <- intersect( feature, feature_names(x) )
        if( feature_are_rows(x) ){
            return(x[feature, , drop=FALSE])
        } else {
            return(x[, feature, drop=FALSE])
        }
    }
})
#' @aliases prune_feature,character,sample_data-method
#' @rdname prune_feature-methods
setMethod("prune_feature", signature("character", "Sample_data"), function(feature, x){
    return(x)
})
#' @aliases prune_feature,character,MyDataSet-method
#' @rdname prune_feature-methods
setMethod("prune_feature", signature("character", "MyDataSet"), function(feature, x){
    # Re-define `feature` as the intersection of feature names for each component AND `feature`
    feature <- intersect(intersect_feature(x), feature)
    # Now prune them all.
    # All MyDataSet objects have an Profile_table slot, no need to test for existence.
    x@Profile_table <- Profile_table(prune_feature(feature, Profile_table(x)), feature_are_rows = TRUE)
    # Test if slot is present. If so, perform the component prune.
    if( !is.null(x@Feature_table) ){
        x@Feature_table <- Feature_table(prune_feature(feature, Feature_table(x)))
    }
    # Force index order after pruning to be the same,
    # according to the same rules as in the constructor, MyDataSet()
    x <- index_reorder(x, index_type="Feature")
    return(x)
})
#' @aliases prune_feature,character,Feature_table-method
#' @rdname prune_feature-methods
setMethod("prune_feature", signature("character", "Feature_table"), function(feature, x){
    if( setequal(feature, feature_names(x)) ){
        return(x)
    } else {
        feature <- intersect( feature, feature_names(x) )
        return( x[feature, , drop=FALSE] )
    }
})
################################################################################
