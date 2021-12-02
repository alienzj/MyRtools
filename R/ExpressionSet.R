#' @importFrom data.table fread
#' @importFrom dplyr %>% intersect select inner_join filter
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom stats setNames
#' @import convert
#'
#'
#' @title Building ExpressionSet Object
#'
#' @description
#' ExpressionSet object comprising not only expressiondata but also metadata is very convenient to perform data analysis.
#'
#' @details 12/2/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param Profile, Numeric matrix; (Required)a Matrix of expression data, whose Row is FeatureID and Column is SampleID.
#' @param Metadata, Data.frame; (Required)a dataframe. of Metadata(1st column must be "SampleID"), containing Group information and also environmental factors(biological factors).
#' @param Feature, Data.frame; the feature of Profile.
#' @param Trim, Logical; filtering Feature and Sample with occurrence(default: FALSE).
#' @param Occ_Feature, Numeric; the cutoff of Feature's trim(default: occurrence=0.2).
#' @param Each, Logical; filtering Features each group or whole data(default: Each("Group")=FALSE).
#' @param Occ_Sample, Numeric; the cutoff of Sample's trim(default: occurrence=0.2).
#'
#' @return
#' an ExpressionSet Object
#'
#' @usage get_ExprSet(profile=Profile, metadata=Metadata, feature=Feature)
#' @examples
#'
#' Profile <- data.table::fread(system.file("extdata", "Species_relative_abundance.tsv", package="MyRtools"))  %>% tibble::column_to_rownames("V1")
#' Metadata <- read.csv(system.file("extdata", "Metadata.csv", package="MyRtools"))
#' Feature <- read.csv(system.file("extdata", "Species_feature.csv", package="MyRtools")) %>% tibble::column_to_rownames("Species")
#'
#' ExprSet <- get_ExprSet(profile=Profile,
#' metadata=Metadata,
#' feature=Feature,
#' trim=TRUE,
#' occ_Feature=0.2,
#' each=TRUE,
#' occ_Sample=0.2)
#'
get_ExprSet <- function(profile=Profile,
                        metadata=Metadata,
                        feature=NULL,
                        trim=FALSE,
                        occ_Feature=0.2,
                        each=FALSE,
                        occ_Sample=0.2){

  # trim feature and sample
  if(trim){
    TrimFun <- function(x, y){
      # All
      feature_occ <- apply(x, 1, function(x){length(x[x!=0])/length(x)}) %>%
        data.frame() %>% stats::setNames("Occ") %>%
        tibble::rownames_to_column("Feature")
      sample_occ <- apply(x, 2, function(x){length(x[x!=0])/length(x)}) %>%
        data.frame() %>% stats::setNames("Occ")

      # Each
      if(each){
        mdat <- y %>% dplyr::select(all_of(c("SampleID", "Group"))) %>%
          dplyr::inner_join(x %>%
                              t() %>% data.frame() %>%
                              tibble::rownames_to_column("SampleID"),
                            by = "SampleID")

        group_name <- unique(mdat$Group)
        group_name_each <- lapply(group_name, function(x){
          mdat_cln <- mdat %>% dplyr::filter(Group%in%x)
          return(mdat_cln$SampleID)
        })

        feature_occ_each <- sapply(1:length(group_name_each), function(i){
          df <- mdat %>% dplyr::filter(SampleID%in%group_name_each[[i]])
          df2 <- df %>% dplyr::select(-"Group") %>%
            tibble::column_to_rownames("SampleID") %>%t()
          ratios <- as.numeric(apply(df2, 1, function(x){length(x[x!=0])/length(x)}))
          return(ratios)
        }) %>% data.frame() %>% stats::setNames(paste0(group_name, "_Occ"))
        rownames(feature_occ_each) <- colnames(mdat)[-c(1:2)]

        feature_occ <- feature_occ_each
      }

      res <- list(feature=feature_occ, sample=sample_occ)
    }

    # Occurrence
    occ_res <- TrimFun(Profile, Metadata)

    # filtering
    feature_KEEP <- apply(occ_res$feature > occ_Feature, 1, all) %>%
      data.frame() %>% stats::setNames("Status") %>%
      dplyr::filter(Status)

    sample_KEEP <- apply(occ_res$sample > occ_Feature, 1, all) %>%
      data.frame() %>% stats::setNames("Status") %>%
      dplyr::filter(Status)

    # run
    profile <- profile[rownames(profile)%in%rownames(feature_KEEP),
                       colnames(profile)%in%rownames(sample_KEEP)]
    message(paste("The remained data is",
                  paste(dim(profile)[1], "Features and",
                        dim(profile)[2], "Samples")))

  }

  intersect_id <- dplyr::intersect(colnames(profile), metadata$SampleID)
  phen <- metadata %>% dplyr::filter(SampleID%in%intersect_id) %>%
    column_to_rownames("SampleID")
  prof <- profile %>% dplyr::select(all_of(rownames(phen)))
  if(!any(rownames(phen) == colnames(prof))){
    stop("Please check the order of SampleID between phen and prof")
  }

  # ExpressionSet
  exprs <- as.matrix(prof)
  adf <- new("AnnotatedDataFrame", data=phen)
  experimentData <- new("MIAME",
                        name="Hua Zou", lab="UCAS",
                        contact="zouhua1@outlook.com",
                        title="Experiment",
                        abstract="Profile",
                        url="www.zouhua.top",
                        other=list(notes="Expression"))
  # Feature input or not
  if(is.null(feature)){
    expressionSet <- new("ExpressionSet",
                         exprs=exprs,
                         phenoData=adf,
                         experimentData=experimentData)
  }else{
    feat <- feature[rownames(feature)%in%rownames(prof), ]
    if(!any(rownames(feat) == rownames(prof))){
      stop("Please check the order of Feature between feat and prof")
    }
    fdf <- new("AnnotatedDataFrame", data=feat)
    expressionSet <- new("ExpressionSet",
                         exprs=exprs,
                         phenoData=adf,
                         featureData=fdf,
                         experimentData=experimentData)
  }

  return(expressionSet)
}
