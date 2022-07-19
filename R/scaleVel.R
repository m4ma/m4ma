#' scaleVel
#' 
#' Scale velocity by time step (tStep)
#' @param v Numeric vector 
#' @param tStep double
#' @return Numeric vector of scaled velocity of same length of v

scaleVel <- function(v, tStep = 0.5) {
  v * tStep  
}