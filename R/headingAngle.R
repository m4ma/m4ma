#' headingAngle
#' 
#' Absolute angular difference between 11 directions with zero cone having angle a1 and a2
#' @param a1 Numeric vector 
#' @param a2 Numeric vector 
#' @return Numeric matrix whose rows and columns are the same length of a2 and a1
headingAngle <- function(a2, a1, angles = c(72.5, 50, 32.5, 20, 10, 0, 350, 
                                            340, 327.5, 310, 287.5)) {
  sapply((angles + a1) %% 360, minAngle, a2 = a2)  
}

