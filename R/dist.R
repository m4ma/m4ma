#' dist
#'
#' Compute distance from p1 to p2, both xy column matrices
#' @param p1 Numeric matrix
#' @param p2 Numeric matrix
#' @return Named numeric vector of length equal to the number of rows in p1
#' @export


dist <- function(p1, p2) {
  sqrt(apply((p1 - p2)^2, 1, sum))
}