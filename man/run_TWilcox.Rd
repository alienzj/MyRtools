% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DA_T_Wilcox_Test.R
\name{run_TWilcox}
\alias{run_TWilcox}
\title{Differential Expression Analysis by T test or Wilcox-rank-sum test}
\usage{
run_Wilcox(dataset=ExpressionSet,
                  trim="none",
                  transform="none",
                  normalize="none",
                  Group_info="Group",
                  Group_name=c("HC", "AA"),
                  Pvalue=0.05,
                  Log2FC=1)
}
\arguments{
\item{trim, }{Character; filter to apply.(default: trim="none").}

\item{transform, }{Character; transformation to apply.(default: tranform="none").}

\item{normalize, }{Character; normalization to apply.(default: normalize="none").}

\item{Group_info, }{Character; design factor(default: "Group").}

\item{Group_name, }{Character; (Required) the group for comparison.}

\item{Pvalue, }{Numeric; significant level(default: 0.05).}

\item{Log2FC, }{Numeric; log2FoldChange(default: 1).}

\item{Expression, }{ExpressionSet; (Required) ExpressionSet object.}
}
\value{
significant difference with enriched directors
}
\description{
The data of features maybe either non-normal distribution or normal distribution(Gauss Distribution).
shapiro.test is used to check the distribution. T-test is for normal distribution, while Wilcox-rank-sum test is for non-normal distribution.
The log2FoldChange is calculate by the geometricmean or median of two groups after a linear scale (Zscore).
}
\details{
12/3/2021 Guangzhou China
}
\examples{

\donttest{
data(ExprSetRawRB)

TWilcox_res <- run_TWilcox(dataset=ExprSetRawRB, Group_info="Group", Group_name=c("HC", "AA"), Pvalue=0.05, Log2FC=1)
TWilcox_res$res
}

}
\author{
Hua Zou
}
