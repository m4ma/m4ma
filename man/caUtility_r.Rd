% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R, R/utility_wrappers.R
\name{caUtility_r}
\alias{caUtility_r}
\alias{caUtility}
\title{Current Angle Utility}
\usage{
caUtility_r(p, angles = c(10, 20, 32.5, 50, 72.5)/90)

caUtility(p, use = "cpp")
}
\arguments{
\item{p}{Numeric vector with subject parameters.}

\item{angles}{Numeric vector with angles of cones.}

\item{use}{Whether the R or C++ implementation is called.}
}
\value{
Numeric vector with current angle utilities for each cell.
}
\description{
Current Angle Utility
}
