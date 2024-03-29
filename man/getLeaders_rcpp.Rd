% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R, R/utility_extra.R,
%   R/utility_extra_wrappers.R
\name{getLeaders_rcpp}
\alias{getLeaders_rcpp}
\alias{getLeaders_r}
\alias{getLeaders}
\title{Leaders}
\usage{
getLeaders_rcpp(
  n,
  p_mat,
  a,
  v,
  P1,
  group,
  centres,
  objects,
  onlyGroup = FALSE,
  preferGroup = TRUE,
  pickBest = FALSE
)

getLeaders_r(
  n,
  state,
  centres,
  objects,
  onlyGroup = FALSE,
  preferGroup = TRUE,
  pickBest = FALSE
)

getLeaders(..., use = "cpp")
}
\arguments{
\item{n}{Integer scalar index of current pedestrian in pedestrian matrix.}

\item{p_mat}{Numeric matrix with shape Nx2 (x and y) indicating the positions of all pedestrians.}

\item{a}{Numeric vector with angles of all pedestrians.}

\item{v}{Numeric vector with velocities of all pedestrians.}

\item{P1}{Numeric matrix with shape Nx2 (x and y) indicating goal positions of current pedestrian.}

\item{group}{Numeric vector with group membership indices of all pedestrians.}

\item{centres}{Numeric matrix with x and y for each cell centre (33x2).}

\item{objects}{List containing a list for each object. An object has
two length-two numeric vectors of x- and y-coordinates.}

\item{onlyGroup}{Logical scalar indicating if leaders must be from the same
group.}

\item{preferGroup}{Logical scalar indicating if leaders are preferred from
the same group.}

\item{pickBest}{Logical salar indicating if all or the best leader should be
returned.}

\item{...}{Arguments passed to the function implementation.}

\item{use}{Whether the R or C++ implementation is called.}
}
\value{
List:
\describe{
\item{dists}{Numeric matrix with distances for each cell (columns) to
each leader (rows).}
\item{leaders}{Numeric matrix with three rows for every leader
(columns). The first row is the cell of the leader, the second the
difference in heading angle between pedestrian \code{n} and leader, the
third the index of the shared group.}
}
}
\description{
Get leaders for pedestrian \code{n}.
}
\details{
Returns \code{NULL} if no other pedestrians are present, if they cannot be seen,
if they are not in front, or if no other pedestrian has the same group as
pedestrian \code{n} and \code{onlyGroup} is \code{TRUE}.
}
