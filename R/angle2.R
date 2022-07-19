#' angle2
#' 
#' Anti-clockwise angle from p1 as origin to p2 (x,y pairs matrices). The angle goes from 0 to 360.
#' @param p1 Numeric matrix
#' @param p2 Numeric matrix
#' @return Named numeric vector of length equal to the number of rows in p1
#' @export


angle2 <- function(p1, p2) {
  #round((atan2(p2[, 2] - p1[, 2], p2[, 1] - p1[, 1]) * 180 / pi), 10) %% 360
  (atan2(p2[, 2] - p1[, 2], p2[, 1] - p1[, 1]) * 180 / pi) %% 360
}
