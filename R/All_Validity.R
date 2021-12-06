######################################################
# https://github.com/joey711/MyDataSet/blob/master/R/
########################################
# Profile_table:
########################################
validProfile_table <- function(object){
    # Both dimensions must have non-zero length.
    if( any(dim(object)==0) ){
        return("\n Profile table data must have non-zero dimensions.")
    }
    # Verify that it is numeric matrix
    if( !is.numeric(object@.Data[, 1]) ){
        text = "\n Non-numeric matrix provided as Profile table.\n"
        text = paste0(text, "Expression is expected to be numeric.")
        return(text)
    }
    return(TRUE)
}
## assign the function as the validity method for the Profile_table class
setValidity("Profile_table", validProfile_table)
########################################
# Sample_data:
########################################
validSample_data <- function(object){
    if( any(dim(object)==0) ){
        return("Sample Data must have non-zero dimensions.")
    }
    return(TRUE)
}
## assign the function as the validity method for the Sample_data class
setValidity("Sample_data", validSample_data)
########################################
# Feature_table:
########################################
validFeature_table <- function(object){
    # Both dimensions must have non-zero length.
    if( any(dim(object)==0) ){
        return("\n Feature Table must have non-zero dimensions.")
    }
    # Verify that it is character matrix
    if( !is.character(object@.Data[, 1]) ){
        text = "\n Non-character matrix provided as Feature Table.\n"
        text = paste0(text, "Feature is expected to be characters.")
        return(text)
    }
    return(TRUE)
}
## assign the function as the validity method for the sample_data class
setValidity("Feature_table", validFeature_table)
########################################
# MyDataSet-class:
########################################
validMyDataSet <- function(object){
    # There must be an Profile_table
    if( is.null(object@Profile_table) ){
        return("\n An Profile_table is required for most analysis / graphics in the MyDataSet-package")
    }
    # intersection of feature-names must have non-zero length
    if( length(intersect_feature(object)) <= 0 ){
        return(paste("\n Component feature names do not match.\n",
                     " Feature indices are critical to analysis.\n Try feature_names()", sep=""))
    }
    # If there is sample data, check that sample-names overlap
    if( !is.null(object@Sample_data) ){
        if( length(intersect(sample_names(object@Sample_data), sample_names(object@Profile_table))) <= 0 ){
            return("\n Component sample names do not match.\n Try sample_names()")
        }
    }
    return(TRUE)
}
setValidity("MyDataSet", validMyDataSet)
########################################
