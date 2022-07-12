#' Walk-beside Utility
#'
#' @param p Numeric vector with subject parameters.
#' @param WB Numeric vector of distances from cells' centers to closest buddy.
#'
#' @return Numeric vector with walk-beside utilities for each cell.
#' @export
#'
wbUtility <- function(p, WB) {
  if (is.null(WB)) {
    return(numeric(33))
  }
  
  # Sum of utilities for buddies
  -apply(sapply(1:dim(WB$dists)[1], function(i) {
    p["bWB"] * WB$buddies["angleDisagree", i] * WB$dists[i, ]^p["aWB"]
  }), 1, sum)
}