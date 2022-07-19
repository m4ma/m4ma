#' aTOd
#' 
#' Compute sine and cosine of an angle.
#' @param a Numeric vector - angles in degrees between 0 and 360
#' @return Numeric Matrix of x and y coordinates (between -1 and 1) that are the signed (i.e., +/-) normalised difference between the xy points that generated the angle a

aTOd <- function(a) {
  cbind(x = cos(a * pi / 180), y = sin(a * pi / 180))
}