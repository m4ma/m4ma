#' minAngle
#' 
#' Shortest absolute angle between a1 and a2
#' @param a1_double scalar vector 
#' @param a2 Numeric vector 
#' @return Numeric vector of length equal to a1
minAngle <- function(a1,a2) {
  pmin(abs(a1 - a2), abs(pmin(a1, a2) + (360 - pmax(a1, a2))))
}