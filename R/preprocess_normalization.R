#' @title Normalize the expression profile data
#'
#' @description
#'
#' Eliminating the bias such as sequencing depth, platform and other technical effects is vital for downstream data analysis.
#'
#' @details 12/7/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @references https://github.com/yiluheihei/microbiomeMarker.
#'
#' @param object, Object; a [`matrix`] or [`assayData-class`] or [`ExpressionSet-class`].
#' @param normalize, Character; normalization to apply, the options inclulde:
#' * "none", return the original data without any transformation.
#' * "TSS": total sum scaling, also referred to as "relative abundance", the
#'     abundances were normalized by dividing the corresponding sample library
#'     size.
#' * "TMM": trimmed mean of m-values. First, a sample is chosen as reference.
#'     The scaling factor is then derived using a weighted trimmed mean over
#'     the differences of the log-transformed gene-count fold-change between
#'     the sample and the reference.
#' * "RLE", relative log expression, RLE uses a pseudo-reference calculated
#'     using the geometric mean of the gene-specific abundances over all
#'     samples. The scaling factors are then calculated as the median of the
#'     gene counts ratios between the samples and the reference.
#' * "CSS": cumulative sum scaling, calculates scaling factors as the
#'     cumulative sum of gene abundances up to a data-derived threshold.
#' * "CLR": centered log-ratio normalization.
#' * "CPM": pre-sample normalization of the sum of the values to 1e+06.
#' * "Zscore": convert the data into unit normal distribution(mean=0; sd=1).
#'
#'
#' @return
#' A object matches the class of argument `object` with the normalized profile.
#'
#' @export
#'
#' @importFrom Biobase exprs assayDataNew
#' @importFrom stats mad median quantile sd
#'
#' @usage run_normalize(object, normalize=c("none", "TSS", "TMM", "RLE", "CSS", "CLR", "CPM", "Zscore"))
#' @examples
#'
run_normalize <- function(object,
                          normalize = c("none", "TSS", "TMM", "RLE", "CSS", "CLR", "CPM", "Zscore")){


  normalize <- match.arg(normalize, c("none", "TSS", "TMM", "RLE", "CSS", "CLR", "CPM", "Zscore"))
  if(inherits(object, "ExpressionSet")){
    prf <- as(Biobase::exprs(object), "matrix")
  }else if(inherits(object, "environment")){
    prf <- as(object$exprs, "matrix")
  }else{
    prf <- object
  }

  if (normalize == "none") {
    abd <- prf
  } else if (normalize == "TSS") {
    abd <- norm_tss(prf)
  } else if (normalize == "TMM")  {
    abd <- norm_tmm(prf)
  }else if (normalize == "RLE") {
    abd <- transform_log10(prf)
  } else if (normalize == "CSS") {
    abd <- transform_log10p(prf)
  } else if (normalize == "CLR") {
    abd <- transform_log10p(prf)
  } else if (normalize == "CPM") {
    abd <- transform_log10p(prf)
  } else if (normalize == "Zscore") {
    abd <- transform_log10p(prf)
  }

  if(inherits(object, "ExpressionSet")){
    object <- get_TransformedExprSet(object, abd)
  }else if(inherits(object, "environment")){
    object <- Biobase::assayDataNew(exprs=abd)
  }else{
    object <- abd
  }

  return(object)
}
################################################################################
#' Total-Sum Scaling (TSS) method
#'
#' TSS simply transforms the feature table into relative abundance by dividing
#' the number of total reads of each sample.
#'
#' @param object, Object; a [`matrix`](Row->Features; Column->Samples).
#'
#' @export
#' @rdname normalize-methods
#' @aliases norm_tss
#'
norm_tss <- function(object){

  size <- colSums(object)
  object_normed <- sweep(object, MARGIN = 2, STATS = size, FUN = "/")

  return(object_normed)
}

################################################################################
# https://github.com/biobakery/Maaslin2/blob/master/R/utility_scripts.R
#
#' TMM (trimmed mean of m-values) normalization
#'
#' TMM calculates the normalization factor using a robust statistics based on
#' the assumption that most features are not differential and should, in
#' average, be equal between the samples. The TMM scaling factor is calculated
#' as the weighted mean of log-ratios between each pair of samples, after
#' excluding the highest count OTUs and OTUs with the largest log-fold change.
#'
#' @param object, Object; a [`matrix`](Row->Features; Column->Samples).
#' @param ref_column column to use as reference
#' @param logratio_trim amount of trim to use on log-ratios
#' @param sum_trim amount of trim to use on the combined absolute levels
#'   ("A" values)
#' @param do_weighting whether to compute the weights or not
#' @param Acutoff cutoff on "A" values to use before trimming
#' @seealso [edgeR::calcNormFactors()]
#' @export
#' @rdname normalize-methods
#' @aliases norm_tmm
#'
norm_tmm <- function(object,
                     ref_column = NULL,
                     logratio_trim = 0.3,
                     sum_trim = 0.05,
                     do_weighting = TRUE,
                     Acutoff = -1e10) {
  nf <- edgeR::calcNormFactors(
    object,
    method = "TMM",
    refcolumn = ref_column,
    logratioTrim = logratio_trim,
    sumTrim = sum_trim,
    doWeighting = do_weighting,
    Acutoff = Acutoff
  )

 nf
}
################################################################################
#' Relative log expression (RLE) normalization
#'
#' RLE assumes most features are not differential and uses the relative
#' abundances to calculate the normalization factor.
#'
#' @param object, Object; a [`matrix`](Row->Features; Column->Samples).
#' @param locfunc a function to compute a location for a sample. By default,
#'   the median is used.
#' @param type method for estimation: either "ratio"or "poscounts" (recommend).
#' @param geo_means default `NULL`, which means the geometric means of the
#'   counts are used. A vector of geometric means from another count matrix can
#'   be provided for a "frozen" size factor calculation.
#' @param control_genes default `NULL`, which means all taxa are used for size
#'   factor estimation, numeric or logical index vector specifying the taxa
#'   used for size factor estimation (e.g. core taxa).
#' @seealso [DESeq2::estimateSizeFactorsForMatrix()]
#' @export
#' @rdname normalize-methods
#' @aliases norm_rle
#'
norm_rle <- function(object,
                     locfunc = stats::median,
                     type = c("poscounts", "ratio"),
                     geo_means = NULL,
                     control_genes = NULL) {


  type <- match.arg(type, c("poscounts", "ratio"))

  # use substitute() to create missing argument
  geo_means <- ifelse(is.null(geo_means), substitute(), geo_means)
  control_genes <- ifelse(is.null(control_genes), substitute(), control_genes)

  nf <- estimateSizeFactorsForMatrix(
    object,
    locfunc = locfunc,
    geoMeans = geo_means,
    controlGenes = control_genes,
    type = type
  )
  object_nf <- set_nf(object, nf)

  object_nf
}
################################################################################
#' Cumulative-Sum Scaling (CSS) method
#'
#' CSS is based on the assumption that the count distributions in each sample
#' are equivalent for low abundant genes up to a certain threshold.  Only the
#' segment of each sample’s count distribution that is relatively invariant
#' across samples is scaled by CSS
#'
#' @param object, Object; a [`matrix`](Row->Features; Column->Samples).
#' @param sl The value to scale.
#' @importFrom phyloseq sample_data<-
#' @importFrom metagenomeSeq newMRexperiment cumNorm cumNormStatFast MRcounts
#' @seealso [metagenomeSeq::calcNormFactors()]
#' @export
#' @rdname normalize-methods
#' @aliases norm_css
#'
norm_css <- function(object, sl = 1000) {

  if (inherits(object, "phyloseq")) {
    object_mgs <- phyloseq2metagenomeSeq(object)
  } else if (inherits(object, "otu_table")) {
    object_mgs <- otu_table2metagenomeSeq(object)
  }

  # cumNormStatFast requires counts of all samples at least have two
  # non zero features. Thus, if there are samples with only one non-zer
  # features, cumNormStat is taken to compute the pth quantile.
  count <- as(otu_table(object), "matrix")
  fun_p <- select_quantile_func(count)
  nf <- metagenomeSeq::calcNormFactors(object_mgs, p = fun_p(object_mgs))
  nf <- unlist(nf) / sl
  object_nf <- set_nf(object, nf)

  object_nf
}
################################################################################
#' CLR (centered log-ratio) normalization
#'
#' In CLR, the log-ratios are computed relative to the geometric mean of all
#' features.
#'
#' @param object, Object; a [`matrix`](Row->Features; Column->Samples).
#'
#' @export
#' @rdname normalize-methods
#' @aliases norm_clr
#'
norm_clr <- function(object) {

  object_normed <- apply(object, 2, trans_clr)

  # do not save the norm_factor, the norm factors are calculated based on the
  # subsequently differential analysis method, e.g. edgeR, DESeq

  return(object_normed)

}

# from joey711/shiny-phyloseq/blob/master/panels/paneldoc/Transform.md
gm_mean <- function(x, na.rm = TRUE) {
  # The geometric mean, with some error-protection bits.
  exp(sum(log(x[x > 0 & !is.na(x)]), na.rm = na.rm) / length(x))
}

trans_clr <- function(x, base = exp(1)) {
  x <- log((x / gm_mean(x)), base)
  x[!is.finite(x) | is.na(x)] <- 0.0
  return(x)
}
################################################################################
#' Normalize the sum of values of each sample to million (counts per million)
#'
#' `norm_cpm`: This normalization method is from the original LEfSe algorithm,
#' recommended when very low values are present (as shown in the LEfSe galaxy).
#'
#' @param object, Object; a [`matrix`](Row->Features; Column->Samples).
#'
#' @export
#' @rdname normalize-methods
#' @aliases norm_cpm
#' @importFrom phyloseq transform_sample_counts
norm_cpm <- function(object) {
  otu <- as(otu_table(object), "matrix") %>%
    as.data.frame()

  # whether the object is summarized
  hie <- check_tax_summarize(object)
  if (hie) {
    features <- row.names(otu)
    features_split <- strsplit(features, "|", fixed = TRUE)
    single_indx <- which(lengths(features_split) < 2)

    ## keep the counts of a sample identical with `normalization`
    ## if we norm the counts in two steps:
    ## 1. calculate scale size: norm_coef = normalization/lib_size;
    ## 2. multiple the scale size value * norm_coef
    ## the counts of a sample colSums(otu) may not equal to the argument
    ## normalization.
    ## e.g. normalization = 1e6, colSums(otu) = 999999
    ## Finally, the kruskal test may be inaccurate,
    ## e.g. https://github.com/yiluheihei/microbiomeMarker/issues/13
    ps_normed <- transform_sample_counts(
      object,
      function(x) x * 1e+06 / sum(x[single_indx])
    )
  } else {
    ps_normed <- transform_sample_counts(
      object,
      function(x) x * 1e+06 / sum(x)
    )
  }

  otu_normed <- data.frame(otu_table(ps_normed))
  otu_normed <- purrr::map_df(
    otu_normed,
    function(x) {
      if (mean(x) && stats::sd(x) / mean(x) < 1e-10) {
        return(round(x * 1e6) / 1e6)
      } else {
        return(x)
      }
    }
  )

  otu_normed <- as.data.frame(otu_normed)
  row.names(otu_normed) <- row.names(otu)
  colnames(otu_normed) <- colnames(otu)
  otu_table(object) <- otu_table(otu_normed, taxa_are_rows = TRUE)

  # do not save the norm_factor, the norm factors are calculated based on the
  # subsequently differential analysis method, e.g. edgeR, DESeq
  object
}
################################################################################
