% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DA_ANCOM.R
\name{run_ANCOM}
\alias{run_ANCOM}
\title{Differential Expression Analysis by Analysis of composition of microbiomes(ANCOM)}
\usage{
run_ANCOM(dataset=ExpressionSet,
                 trim="none",
                 GroupVar="Group",
                 AdjVar=NULL,
                 RandVar=NULL,
                 Pvalue=0.05,
                 Wvalue=0.7)
}
\arguments{
\item{trim, }{Character; filter to apply.(default: trim="none").}

\item{GroupVar, }{Character; design factor(default: "Group").}

\item{AdjVar, }{Character; the adjusted variables.}

\item{RandVar, }{Character; random effects for longitudinal analysis or repeated measure("~ 1 | studyid").}

\item{Pvalue, }{Numeric; significant level(default: 0.05).}

\item{Wvalue, }{Numeric; W statistic for clarify the significant features(default: 0.7).}

\item{Expression, }{ExpressionSet; (Required) ExpressionSet object.}
}
\value{
a list of results:
  significant features pass the threshold of W statistics
  A volcano plot of significant features
}
\description{
ANCOM accounts for the underlying structure in the data and can be used for comparing
the composition of microbiomes in two or more populations.
ANCOM makes no distributional assumptions and can be implemented in
a linear model framework to adjust for covariates as well as model longitudinal data.
ANCOM also scales well to compare samples involving thousands of taxa.
}
\details{
12/3/2021 Guangzhou China
}
\examples{

\donttest{
data(ExprSetRawCount)

ANCOM_res <- run_ANCOM(dataset=ExprSetRawCount, GroupVar="Group", AdjVar="Gender", RandVar=NULL, Pvalue=0.05, Wvalue=0.7)
ANCOM_res$res
}

}
\references{
Mandal et al. "Analysis of composition of microbiomes: a novel
method for studying microbial composition", Microbial Ecology in Health
& Disease, (2015), 26.
}
\author{
Hua Zou
}
