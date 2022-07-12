#' Preferred Speed Utility
#'
#' @param p Numeric vector with subject parameters.
#' @param v Numeric scalar indicating the current speed.
#' @param d Numeric scalar indicating the distance to next goal.
#'
#' @return Numeric vector with preferred speed utilities for each cell.
#' @export
#'
psUtility <- function(p, v, d) {
  sPref <- pmin(p["sPref"], d * p["sPref"] / p["sSlow"]) # step 1
  -p["bPS"] * c(rep(abs(v * 1.5 - sPref)^p["aPS"], 11), 
                rep(abs(v - sPref)^p["aPS"],11),
                rep(abs(v/2 - sPref)^p["aPS"], 11))  
}