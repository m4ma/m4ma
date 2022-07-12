#' Goal Angle Utility
#'
#' @param p Numeric vector with subject parameters.
#' @param GA Numeric vector with angles to next goal.
#'
#' @return Numeric vector with goal angle utilities for each cell.
#' @export
#'
gaUtility <- function(p, GA) {
  -rep(p["bGA"] * GA^p["aGA"], times = 3)
}

