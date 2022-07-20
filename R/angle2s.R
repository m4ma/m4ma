#'
#' Compute Shortest angle anti-clockwise from p1 as origin to p2 (> -180 to 180)
#' @param p1 Numeric matrix
#' @param p2 Numeric matrix
#' @return Numeric vector of length equal to the number of rows in p1 
#' @export
angle2s <- function(p1, p2) {
  round((atan2(p2[, 2] - p1[, 2], p2[, 1] - p1[, 1]) * 180 / pi), 10)
}
