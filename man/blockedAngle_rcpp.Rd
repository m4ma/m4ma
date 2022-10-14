% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R, R/utility_extra.R,
%   R/utility_extra_wrappers.R
\name{blockedAngle_rcpp}
\alias{blockedAngle_rcpp}
\alias{blockedAngle_r}
\alias{blockedAngle}
\title{Blocked Angle}
\usage{
blockedAngle_rcpp(n, state, p_pred, objects)

blockedAngle_r(n, state, p_pred, objects)

blockedAngle(..., use = "cpp")
}
\arguments{
\item{n}{Integer scalar index of current pedestrian in pedestrian matrix.}

\item{state}{List with state data.}

\item{p_pred}{Numeric matrix with shape Nx2 (x and y) as predicted positions
of all pedestrians.}

\item{objects}{List containing a list for each object. An object has
two length-two numeric vectors of x- and y-coordinates.}

\item{...}{Arguments passed to the function implementation.}

\item{use}{Whether the R or C++ implementation is called.}
}
\value{
Numeric matrix with rows for each other pedestrian and columns for
each cell (Nx33) or NULL.
}
\description{
Compute the distance from each cell to closest pedestrians.
}