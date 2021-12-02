#' @importFrom stats mad median quantile sd
#'
#'
#' @title Normalization
#'
#' @description
#' Multiple variables may have different, leading to unequal contribution on the data analysis. Normalization could figure out this drawback when we do some predictive analysis.
#'
#' @details 12/1/2021 Guangzhou China
#' @author  Hua Zou
#'
#'
#' @param Vectors, a vector contains numeric.
#' @param type, the default method is Zscore. There are some other approaches, such as "MAD": Median Absolute Deviation;
#' "Robust"; "Unit"; "Zscore"; "Median"; "Min_Max" scale normalization.
#'
#' @return
#' a normalized vector
#'
#' @usage NormalizeFun(Vectors, type="Zscore")
#' @examples
#'
#' V <- c(1:20, 30:40)
#' M <- matrix(rnorm(24, mean = 10, sd=10), nrow=4, ncol=6, byrow=TRUE, dimnames = list(paste0("R", 1:4), paste0("C", 1:6)))
#' V_norm <- NormalizeFun(V, type="Zscore")
#' M_norm <- t(apply(M, 1, NormalizeFun, "Zscore"))

NormalizeFun <- function(Vectors, type="Zscore"){


  # Median Absolute Deviation normalization
  MAD_norm <- function(features){
    # x for features X = (x1, x2, ..., xn)
    value <- as.numeric(features)
    d_mad <- stats::mad(value, na.rm = FALSE)
    x_scale <- (value - stats::median(value, na.rm = FALSE))/d_mad

    return(x_scale)
  }

  # Robust scale normalization
  Robust_norm <- function(features){
    # x for features X = (x1, x2, ..., xn)
    value <- as.numeric(features)
    q_value <- as.numeric(stats::quantile(value, na.rm = FALSE))
    remain_value <- value[value > q_value[2] & value < q_value[4]]

    mean_value <- mean(remain_value, na.rm = FALSE)
    sd_value <- stats::sd(remain_value, na.rm = FALSE)

    x_scale <- (value - mean_value)/sd_value

    return(x_scale)
  }

  # Unit scale normalization
  Unit_norm <- function(samples){
    # v for samples v = (v1, v2, ..., vn)
    value <- as.numeric(samples)
    x_scale <- value / sqrt(sum(value^2, na.rm = FALSE))

    return(x_scale)
  }

  # Z-scale normalization
  Zscore_norm <- function(features){
    # x for features X = (x1, x2, ..., xn)
    value <- as.numeric(features)
    mean_value <- mean(value, na.rm = FALSE)
    sd_value <- stats::sd(value, na.rm = FALSE)

    x_scale <- (value - mean_value)/sd_value

    return(x_scale)
  }

  # Median-scale normalization
  Median_norm <- function(features){
    # x for features X = (x1, x2, ..., xn)
    value <- as.numeric(features)
    median_value <- stats::median(value, na.rm = FALSE)
    IQR_value <- stats::IQR(value, na.rm = FALSE)

    x_scale <- (value - median_value)/IQR_value

    return(x_scale)
  }

  # Min-Max  normalization
  Min_Max_norm <- function(features){
    # x for features X = (x1, x2, ..., xn)
    value <- as.numeric(features)
    min_value <- min(value, na.rm = FALSE)
    max_value <- max(value, na.rm = FALSE)

    x_scale <- (value - min_value)/(max_value - min_value)

    return(x_scale)
  }

  MethodSet <- c("MAD", "Robust", "Unit", "Zscore", "Median", "Min_Max")
  if(is.element(type, MethodSet)){
    if(type == MethodSet[1]){
      res <- MAD_norm(Vectors)
    }else if(type == MethodSet[2]){
      res <- Robust_norm(Vectors)
    }else if(type == MethodSet[3]){
      res <- Unit_norm(Vectors)
    }else if(type == MethodSet[4]){
      res <- Zscore_norm(Vectors)
    }else if(type == MethodSet[5]){
      res <- Median_norm(Vectors)
    }else if(type == MethodSet[6]){
      res <- Min_Max_norm(Vectors)
    }
  }else{
    stop("No methods have been chosen, please check your input")
  }

  return(res)
}
