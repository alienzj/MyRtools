#' @title Building ExpressionSet Object
#'
#' @description
#' ExpressionSet object comprising not only expressiondata but also metadata is very convenient to perform data analysis.
#'
#' @details 12/2/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param profile, Numeric matrix; (Required)a Matrix of expression data, whose Row is FeatureID and Column is SampleID.
#' @param metadata, Data.frame; (Required)a dataframe. of Metadata(1st column must be "SampleID"), containing Group information and also environmental factors(biological factors).
#' @param feature, Data.frame; the feature of Profile(default: Feature=NULL).
#' @param name, Character; name of experimentData, inheriting from ExpressionSet(default: name="").
#' @param lab, Character; lab of experimentData(default: lab="").
#' @param contact, Character; contact of experimentData(default: contact="").
#' @param title, Character; title of experimentData(default: title="").
#' @param abstract, Character; abstract of experimentData(default: abstract="").
#' @param url, Character; url of experimentData(default: url="").
#' @param notes, Character; notes of experimentData(default: notes="").
#'
#' @return
#' an ExpressionSet Object
#'
#' @export
#'
#' @importFrom data.table fread
#' @importFrom dplyr %>% intersect select inner_join filter all_of
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom stats setNames
#' @import Biobase
#'
#' @usage get_ExprSet(profile=Profile,
#'                   metadata=Metadata,
#'                   feature=Feature,
#'                   name="",
#'                   lab="",
#'                   contact="",
#'                   title="",
#'                   abstract="",
#'                   url="",
#'                   notes="")
#' @examples
#'
#' \donttest{
#' library(dplyr)
#' Profile <- data.table::fread(system.file("extdata", "Species_relative_abundance.tsv", package="MyRtools"))  %>% tibble::column_to_rownames("V1")
#' Metadata <- read.csv(system.file("extdata", "Metadata.csv", package="MyRtools"))
#' Feature <- read.csv(system.file("extdata", "Species_feature.csv", package="MyRtools")) %>% tibble::column_to_rownames("Species")
#'
#' ExprSet <- get_ExprSet(profile=Profile,
#'                        metadata=Metadata,
#'                        feature=Feature,
#'                        name="Hua Zou",
#'                        lab="UCAS",
#'                        contact="zouhua1@outlook.com",
#'                        title="Experiment",
#'                        abstract="Profile",
#'                        url="www.zouhua.top",
#'                        notes="Expression")
#' }
#'
get_ExprSet <- function(profile=Profile,
                        metadata=Metadata,
                        feature=NULL,
                        name="",
                        lab="",
                        contact="",
                        title="",
                        abstract="",
                        url="",
                        notes=""){

  shared_samples <- dplyr::intersect(colnames(profile), metadata$SampleID)
  if(any(length(shared_samples) == 0, length(shared_samples) == 1)){
    stop("There are no common samples or just one sample between profile and metadata please check input data\n")
  }
  phen <- metadata %>% dplyr::filter(SampleID%in%shared_samples) %>%
    column_to_rownames("SampleID")
  prof <- profile %>% dplyr::select(dplyr::all_of(rownames(phen)))
  if(!any(rownames(phen) == colnames(prof))){
    stop("Please check the order of SampleID between phen and prof")
  }

  # ExpressionSet
  exprs <- as.matrix(prof)
  adf <- new("AnnotatedDataFrame", data=phen)
  experimentData <- new("MIAME",
                        name=name, lab=lab,
                        contact=contact,
                        title=title,
                        abstract=abstract,
                        url=url,
                        other=list(notes=notes))
  # Feature input or not
  if(is.null(feature)){
    expressionSet <- new("ExpressionSet",
                         exprs=exprs,
                         phenoData=adf,
                         experimentData=experimentData)
  }else{
    shared_feature <- dplyr::intersect(rownames(prof), rownames(feature))
    if(any(length(shared_feature) == 0, length(shared_feature) == 1)){
      stop("There are no common features or just one feature between profile and metadata please check input data\n")
    }
    fdata <- feature[rownames(feature)%in%shared_feature, ]
    prof2 <- prof %>% t() %>% data.frame() %>%
      dplyr::select(dplyr::all_of(rownames(fdata))) %>%
      t() %>% data.frame()
    if(!any(rownames(fdata) == rownames(prof2))){
      stop("Please check the order of Feature between feat and prof")
    }

    fdf <- new("AnnotatedDataFrame", data=fdata)
    exprs <- as.matrix(prof2)
    expressionSet <- new("ExpressionSet",
                         exprs=exprs,
                         phenoData=adf,
                         featureData=fdf,
                         experimentData=experimentData)
  }

  return(expressionSet)
}

#' @title Building ExpressionSet Object with preprocess
#'
#' @description
#' precessing the ExpressionSet object by filtering, transforming or normalization
#'
#' @details 12/7/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param dataset, ExpressionSet; (Required) an Raw ExpressionSet object by `get_ExprSet`.
#' @param trim_cutoff, Numeric; the threshold for filtering (default: Cutoff=0.2).
#' @param trim, Character; filter to apply.(default: trim="none").
#' @param tranform, Character; transformation to apply.(default: tranform="none").
#' @param normalize, Character; normalization to apply.(default: normalize="none").
#' @param inpute, Character; inputation to apply.(default: inpute="none").
#'
#' @return
#' an preprocessed ExpressionSet Object
#'
#' @export
#'
#' @importFrom data.table fread
#' @importFrom dplyr %>% intersect select inner_join filter all_of
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom stats setNames
#' @import Biobase
#'
#' @usage get_processedExprSet(dataset=ExpressionSet,
#'         trim_cutoff=0.2,
#'         trim=c("none", "both", "feature", "sample", "Group"),
#'         transform=c("none", "log2", "log2p", "log10", "log10p"),
#'         normalize=c("none", "TSS", "TMM", "RLE", "CLR", "Zscore", "Median", "MAD", "Robust", "Unit", "Min_Max"),
#'         impute=c("none", "deletion", "GlobalStandard", "KNN"))
#' @examples
#'
#' \donttest{
#' library(dplyr)
#' data("ExprSetRawRB")
#' ExprSet <- get_processedExprSet(dataset=ExprSetRawCount,
#' trim_cutoff=0.2,
#' trim="Group",
#' transform="log2",
#' normalize=NULL,
#' impute=NULL)
#' }
#'
get_processedExprSet <- function(dataset=ExprSetRawCount,
                                 trim_cutoff=0.2,
                                 trim="none",
                                 transform="none",
                                 normalize=NULL,
                                 impute=NULL){

  # trim features and samples
  if(!is.null(trim_cutoff)){
    trim <- match.arg(trim, c("none", "both", "feature", "sample", "Group"))
    filterObject <- run_trim(dataset, trim_cutoff, trim)
    tempObject <- get_TrimedExprSet(filterObject)
  }else{
    tempObject <- dataset
  }

  # transform expression matrix
  if(!is.null(transform)){
    transform <- match.arg(transform, c("none", "log2", "log2p", "log10", "log10p"))
    transformObject <- run_transform(tempObject, transform)
    tempObject <- get_TransformedExprSet(transformObject, exprs(transformObject))
  }else if(any(!is.null(trim), !is.null(transform))){
    tempObject <- tempObject
  }else if(all(!is.null(trim), !is.null(transform))){
    tempObject <- dataset
  }

  # Normalize expression matrix
  if(!is.null(normalize)){
    normalize <- match.arg(normalize, c("none", "TSS", "TMM", "RLE", "CLR", "Zscore", "Median", "MAD", "Robust", "Unit", "Min_Max"))
    normalizeObject <- run_normalize(tempObject, normalize)
    tempObject <- get_NormalizedExprSet(transformObject, exprs(normalizeObject))
  }else if(any(!is.null(trim), !is.null(transform), !is.null(normalize))){
    tempObject <- tempObject
  }else if(all(!is.null(trim), !is.null(transform), !is.null(normalize))){
    tempObject <- dataset
  }

  # impute expression matrix for missing value(NAs) or Zero values
  if(!is.null(impute)){
    impute <- match.arg(impute, c("none", "deletion", "GlobalStandard", "KNN"))
    imputeObject <- run_impute(tempObject, impute)
    tempObject <- get_imputedExprSet(imputeObject, exprs(imputeObject))
  }else if(any(!is.null(trim), !is.null(transform),
               !is.null(normalize), !is.null(impute))){
    tempObject <- tempObject
  }else if(all(!is.null(trim), !is.null(transform),
               !is.null(normalize), !is.null(impute))){
    tempObject <- dataset
  }

  return(tempObject)
}
