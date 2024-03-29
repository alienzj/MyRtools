% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AA_SimpleCorrelation.R
\name{AA_SimpleCC}
\alias{AA_SimpleCC}
\title{Association Analysis by Spearman correlation coefficient}
\usage{
AA_SimpleCC(dataset=ExpressionSet, Group_info="Group", Group_name=NULL, Object="Feature")
}
\arguments{
\item{Group_info, }{Character; design factor(default: "Group").}

\item{Group_name, }{Character; the group for filtering(default: NULL).}

\item{Object, }{Character; (Required) the variables for correlation(Feature;Meta;Mix).}

\item{Expression, }{ExpressionSet; (Required) ExpressionSet object.}
}
\value{
a list of results
  Spearman Correlation Coefficient matrix
  Spearman Correlation Coefficient table
}
\description{
The correlation among features is calculated by Spearman.
}
\details{
12/3/2021 Guangzhou China
}
\examples{

\donttest{
data(ExprSetRawRB)

SimpleCC_res <- AA_SimpleCC(dataset=ExprSetRawRB, Group_info="Group", Group_name=NULL, Object="Feature")

}

}
\author{
Hua Zou
}
