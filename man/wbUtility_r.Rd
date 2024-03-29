% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R, R/utility_wrappers.R
\name{wbUtility_r}
\alias{wbUtility_r}
\alias{wbUtility}
\title{Walk-beside Utility}
\usage{
wbUtility_r(p, WB)

wbUtility(p, WB, use = "cpp")
}
\arguments{
\item{p}{Numeric vector with subject parameters.}

\item{WB}{Numeric vector of distances from cells' centers to closest buddy.}

\item{use}{Whether the R or C++ implementation is called.}
}
\value{
Numeric vector with walk-beside utilities for each cell.
}
\description{
Walk-beside Utility
}
