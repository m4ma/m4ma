#' Follow-leader Utility
#'
#' @param p Numeric vector with subject parameters.
#' @param FL List of two numeric matrices:
#' \describe{
#'   \item{leaders}{Matrix with columns per leader and rows of their normalized angle disagreement and in-group status.}
#'   \item{dists}{Matrix with rows per leader and columns per cell with distances from each cell to chosen cell.}
#' }
#' @return Numeric vector with follow-leader utilities for each cell.
#' @export
#' 
flUtility <- function(p, FL) {
  if (is.null(FL)) {
    return(numeric(33))
  }
  
  # Ingroup and same direction weighted b for each leader
  b <- (p["bFL"] + p["dFL"] * FL$leaders["inGroup", ]) * 
    FL$leaders["angleDisagree", ]
  
  # Sum of utilities for leaders
  -apply(sapply(1:length(b), function(i) {
    b[i] * FL$dists[i, ]^p["aFL"]
  }), 1, sum)
}