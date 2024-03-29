% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess_ExpressionSet.R
\name{get_processedExprSet}
\alias{get_processedExprSet}
\title{Building ExpressionSet Object with preprocess}
\usage{
get_processedExprSet(dataset=ExpressionSet,
        trim_cutoff=0.2,
        trim=c("none", "both", "feature", "sample", "Group"),
        transform=c("none", "log2", "log2p", "log10", "log10p"),
        normalize=c("none", "TSS", "TMM", "RLE", "CLR", "Zscore", "Median", "MAD", "Robust", "Unit", "Min_Max"),
        impute=c("none", "deletion", "GlobalStandard", "KNN"))
}
\arguments{
\item{dataset, }{ExpressionSet; (Required) an Raw ExpressionSet object by `get_ExprSet`.}

\item{trim_cutoff, }{Numeric; the threshold for filtering (default: Cutoff=0.2).}

\item{trim, }{Character; filter to apply.(default: trim="none").}

\item{normalize, }{Character; normalization to apply.(default: normalize="none").}

\item{tranform, }{Character; transformation to apply.(default: tranform="none").}

\item{inpute, }{Character; inputation to apply.(default: inpute="none").}
}
\value{
an preprocessed ExpressionSet Object
}
\description{
precessing the ExpressionSet object by filtering, transforming or normalization
}
\details{
12/7/2021 Guangzhou China
}
\examples{

\donttest{
library(dplyr)
data("ExprSetRawRB")
ExprSet <- get_processedExprSet(dataset=ExprSetRawCount,
trim_cutoff=0.2,
trim="Group",
transform="log2",
normalize=NULL,
impute=NULL)
}

}
\author{
Hua Zou
}
