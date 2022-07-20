#' cones number
#' @param k in 1..33
#' @return Numeric vector of cone numbers of length equal to 33
#' @export
coneNum <- function(k) {
  1 + ((k-1) %% 11)
}

#' ring number
#' @param k in 1..33
#' @return Numeric vector of ring numbers of length equal to 33
#' @export

ringNum <- function(k) {
  1 + (k-1) %/% 11
}
