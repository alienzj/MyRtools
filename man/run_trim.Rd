% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess_trim.R
\name{run_trim}
\alias{run_trim}
\title{Trim samples or features in `assayData` by Occurrences}
\usage{
run_trim(object, object=0.2, trim="Both")
}
\arguments{
\item{object, }{Object; a [`matrix`] or [`assayData-class`] or [`ExpressionSet-class`].}

\item{cutoff, }{Numeric; the threshold for triming (default: Cutoff=0.2).}

\item{trim, }{Character; transformation to apply, the options inclulde:
* "none", return the original data without any trim.
* "both", trim the features and samples of profile more than cutoff.
* "feature", trim the features of profile more than cutoff.
* "sample", trim the samples of profile more than cutoff.
* "Group", trim the features per Group, the samples of the profile more than cutoff, respectively.}
}
\value{
A trimed `object` with profile more than cutoff.
}
\description{
trim samples or features in `assayData` by Occurrences,
which means the samples or features will be discard if they could not pass the cutoff.
}
\details{
12/6/2021 Guangzhou China
}
\examples{
\donttest{
   data("ExprSet_species")
   run_trim(ExprSet_species, object=0.2, trim="both")
}

}
\author{
Hua Zou
}
