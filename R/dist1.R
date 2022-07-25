#' dist1
#'
#' Compute distance from p1 to p2
#' @param p1 Numeric matrix of a single point (i.e., 1 row, 2 xy columns)
#' @param p2 Numeric matrix of multiple xy points
#' @return Named numeric vector of length equal to the number of rows in p2
#' @export

dist1 <- function(p1, p2) {
  sqrt(apply((t(p2) - as.vector(p1))^2, 2, sum))
}