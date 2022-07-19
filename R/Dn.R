#' Dn
#' 
#' Anti-clockwise angle to destination for all pedestrians. The angle goes from 0 to 360
#' @param p_n Numeric matrix of x and y coordinates
#' @param P_n Numeric matrix of x and y coordinates
#' @return Named numeric vector of length equal to the number of rows in p1
Dn <- function(p_n, P_n) {
  out <- numeric(dim(p_n)[1])
  for (i in 1:dim(p_n)[1]) {
    out[i] <- angle2(p_n[i, , drop = FALSE], P_n[i, , drop = FALSE])
  }
  out
}