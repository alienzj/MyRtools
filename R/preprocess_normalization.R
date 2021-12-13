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
#' * "CLR": centered log-ratio normalization.
#' * "Zscore": convert the data into unit normal distribution(mean=0; sd=1).
#' * "Median": Median scale normalization.
#' * "MAD": Median Absolute Deviation.
#' * "Robust": Robust scale normalization.
#' * "Unit": Unit scale normalization.
#' * "Min_Max": Median-scale scale normalization.
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
#' @usage run_normalize(object, normalize=c("none", "TSS", "TMM", "RLE", CLR", "Zscore", "Median", "MAD", "Robust", "Unit", "Min_Max"))
#' @examples
#'
run_normalize <- function(object,
                          normalize = c("none", "TSS", "TMM", "RLE", "CLR", "Zscore", "Median", "MAD", "Robust", "Unit", "Min_Max")){

  data("ExprSetRawCount")
  object=ExprSetRawCount
  normalize = "TSS"

  normalize <- match.arg(normalize, c("none", "TSS", "TMM", "RLE", "CLR", "Zscore", "Median", "MAD", "Robust", "Unit", "Min_Max"))
  if(inherits(object, "ExpressionSet")){
    prf <- as(Biobase::exprs(object), "matrix")
  }else if(inherits(object, "environment")){
    prf <- as(object$exprs, "matrix")
  }else{
    prf <- object
  }

  if (normalize == "none"){
    abd <- prf
  } else if (normalize == "TSS"){
    abd <- norm_tss(prf)
  } else if (normalize == "TMM"){
    abd <- norm_tmm(prf)
  }else if (normalize == "RLE"){
    abd <- norm_rle(prf)
  }else if (normalize == "CLR"){
    abd <- norm_clr(prf)
  }else if (normalize == "Zscore"){
    abd <- norm_zscore(prf)
  }else if (normalize == "Median"){
    abd <- norm_median(prf)
  }else if (normalize == "MAD"){
    abd <- norm_mad(prf)
  }else if (normalize == "Robust"){
    abd <- norm_robust(prf)
  }else if (normalize == "Unit"){
    abd <- norm_unit(prf)
  }else if (normalize == "Min_Max"){
    abd <- norm_minmax(prf)
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

  return(t(t(object)/nf))
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

  nf <- DESeq2::estimateSizeFactorsForMatrix(
    object,
    locfunc = locfunc,
    geoMeans = geo_means,
    controlGenes = control_genes,
    type = type
  )

  return(t(t(object)/nf))
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

  rownames(object_normed) <- rownames(object)
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
#' Z-scale normalization
#'
#' @param object, Object; a [`matrix`](Row->Features; Column->Samples).
#'
#' @export
#' @rdname normalize-methods
#' @aliases norm_zscore
#'
norm_zscore <- function(object){

  trans_zscore <- function(x){
    value <- as.numeric(x)
    mean_value <- mean(value, na.rm = FALSE)
    sd_value <- stats::sd(value, na.rm = FALSE)

    x_scale <- (value - mean_value)/sd_value

    return(x_scale)
  }
  object_normed <- apply(object, 2, trans_zscore)
  rownames(object_normed) <- rownames(object)

  return(object_normed)
}
################################################################################
#' Median-scale normalization
#'
#' @param object, Object; a [`matrix`](Row->Features; Column->Samples).
#'
#' @export
#' @rdname normalize-methods
#' @aliases norm_median
#'
norm_median <- function(object){

  trans_median <- function(x){
    value <- as.numeric(x)
    median_value <- stats::median(value, na.rm = FALSE)
    IQR_value <- stats::IQR(value, na.rm = FALSE)

    x_scale <- (value - median_value)/IQR_value

    return(x_scale)
  }
  object_normed <- apply(object, 2, trans_median)
  rownames(object_normed) <- rownames(object)

  return(object_normed)
}
################################################################################
#' Median Absolute Deviation normalization
#'
#' @param object, Object; a [`matrix`](Row->Features; Column->Samples).
#'
#' @export
#' @rdname normalize-methods
#' @aliases norm_mad
#'
norm_mad <- function(object){

  trans_mad <- function(x){
    value <- as.numeric(x)
    d_mad <- stats::mad(value, na.rm = FALSE)
    x_scale <- (value - stats::median(value, na.rm = FALSE))/d_mad

    x_scale <- (value - median_value)/IQR_value

    return(x_scale)
  }
  object_normed <- apply(object, 2, trans_median)
  rownames(object_normed) <- rownames(object)
  return(object_normed)
}
################################################################################
#' Robust scale normalization
#'
#' @param object, Object; a [`matrix`](Row->Features; Column->Samples).
#'
#' @export
#' @rdname normalize-methods
#' @aliases norm_robust
#'
norm_robust <- function(object){

  trans_robust <- function(x){

    value <- as.numeric(x)
    q_value <- as.numeric(stats::quantile(value, na.rm = FALSE))
    remain_value <- value[value > q_value[2] & value < q_value[4]]

    mean_value <- mean(remain_value, na.rm = FALSE)
    sd_value <- stats::sd(remain_value, na.rm = FALSE)

    x_scale <- (value - mean_value)/sd_value

    return(x_scale)
  }
  object_normed <- apply(object, 2, trans_robust)
  rownames(object_normed) <- rownames(object)

  return(object_normed)
}
################################################################################
#' Unit scale normalization
#'
#' @param object, Object; a [`matrix`](Row->Features; Column->Samples).
#'
#' @export
#' @rdname normalize-methods
#' @aliases norm_unit
#'
norm_unit <- function(object){

  trans_unit <- function(x){
    value <- as.numeric(x)
    x_scale <- value / sqrt(sum(value^2, na.rm = FALSE))

    return(x_scale)
  }
  object_normed <- apply(object, 2, trans_unit)
  rownames(object_normed) <- rownames(object)
  return(object_normed)
}
################################################################################
#' Min-Max  normalization
#'
#' @param object, Object; a [`matrix`](Row->Features; Column->Samples).
#'
#' @export
#' @rdname normalize-methods
#' @aliases norm_minmax
#'
norm_minmax <- function(object){

  trans_minmax <- function(x){
    value <- as.numeric(x)
    min_value <- min(value, na.rm = FALSE)
    max_value <- max(value, na.rm = FALSE)

    x_scale <- (value - min_value)/(max_value - min_value)

    return(x_scale)
  }
  object_normed <- apply(object, 2, trans_minmax)
  rownames(object_normed) <- rownames(object)
  return(object_normed)
}
################################################################################
