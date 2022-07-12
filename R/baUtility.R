#' Blocked Angle Utility
#'
#' @param p Numeric vector with subject parameters.
#' @param BA Numeric vector of distances from each cell to closest pedestrian.
#'
#' @return Numeric vector with blocked angle utilities for each cell.
#' @export
#'
baUtility <- function(p, BA) {
  if (is.null(BA)) {
    return(rep(0,33))
  }
  cells <- rep(0, 33)
  cells[as.numeric(names(BA))] <- p["bBA"] / (pmax(BA, 0)^p["aBA"])
  return(-cells)
}