% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Global_Ordination.R
\name{run_PCoA}
\alias{run_PCoA}
\title{Dimensionality Reduction Analysis: Principal Coordinate Analysis(PCoA)}
\usage{
run_PCoA(dataset=ExpressionSet, Method="bray", Group_info="Group")
}
\arguments{
\item{Method, }{Character; method to calculate the distance(default: bray).}

\item{Group_info, }{Character; the group for plot(default: "Group").}

\item{Expression, }{ExpressionSet; (Required) ExpressionSet object.}
}
\value{
a list object:
  PCoA score
  Result of PERMANOVA
}
\description{
PCoA uses the distance among samples, which calculates through the multiple variables.
}
\details{
12/2/2021 Guangzhou China
}
\examples{

\donttest{
data(ExprSetRawRB)

PCoA_res <- run_PCoA(dataset=ExprSetRawRB, Method="bray", Group_info="Group")
PCoA_res$PCoA
}

}
\author{
Hua Zou
}
