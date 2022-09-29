
#' Blocked Angle Utility
#'
#' @param p Numeric vector with subject parameters.
#' @param BA Numeric vector of distances from each cell to closest pedestrian.
#'
#' @return Numeric vector with blocked angle utilities for each cell.
#' @export
#'
baUtility_r <- function(p, BA) {
  if (is.null(BA)) {
    return(rep(0,33))
  }
  cells <- rep(0, 33)
  cells[as.numeric(names(BA))] <- p["bBA"] / (pmax(BA, 0)^p["aBA"])
  return(-cells)
}


#' Current Angle Utility
#'
#' @param p Numeric vector with subject parameters.
#' @param angles Numeric vector with angles of cones.
#'
#' @return Numeric vector with current angle utilities for each cell.
#' @export
#'
caUtility_r <- function(p, angles = c(10, 20, 32.5, 50, 72.5) / 90) {
  ap <- angles^p["aCA"]
  -rep(c(rep(p["bCA"] * p["bCAlr"], 5), 1, rep(p["bCA"] / p["bCAlr"], 5)) * 
         c(ap[5:1], 0, ap), times = 3)
}


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
flUtility_r <- function(p, FL) {
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


#' Goal Angle Utility
#'
#' @param p Numeric vector with subject parameters.
#' @param GA Numeric vector with angles to next goal.
#'
#' @return Numeric vector with goal angle utilities for each cell.
#' @export
#'
gaUtility_r <- function(p, GA) {
  -rep(p["bGA"] * GA^p["aGA"], times = 3)
}


#' Interpersonal Distance Utility
#'
#' @param p Numeric vector with subject parameters.
#' @param n Integer scalar indexing the subject in the state.
#' @param ID Numeric matrix of predicted distances from the subject to other pedestrians in the front.
#' @param ok Logical matrix indicating if cells are blocked.
#' @param group Integer vector with group indices for each pedestrian.
#'
#' @return Numeric vector with interpersonal distance utilities for each cell.
#' @export
#'
idUtility_r <- function(p, n, ID, ok, group) {
  # None in front return -Inf for cells blocked by objects
  if (is.null(ID) | !any(ok)) {  
    return(as.vector(ifelse(ok, 0, -Inf)))
  }
  
  # Group dependent b, bigger for outgroup by dID
  namesInGroup <- names(group[-n][group[-n] == group[n]])
  bID <- ifelse(dimnames(ID)[[1]] %in% namesInGroup, p["bID"], 
                p["bID"] + p["dID"])
  
  # Object or ped clash
  ok <- ok & apply(ID, 2, function(x) {
    all(x > 0)
  })
  out <- ifelse(ok, 0, -Inf)
  
  # Repulsion
  if (p["bID"] != 0) {
    out[ok] <- -apply(bID / (ID[, ok, drop = FALSE]^p["aID"]), 2, sum)
  }
  
  as.vector(out)
}


#' Preferred Speed Utility
#'
#' @param p Numeric vector with subject parameters.
#' @param v Numeric scalar indicating the current speed.
#' @param d Numeric scalar indicating the distance to next goal.
#'
#' @return Numeric vector with preferred speed utilities for each cell.
#' @export
#'
psUtility_r <- function(p, v, d) {
  sPref <- pmin(p["sPref"], d * p["sPref"] / p["sSlow"]) # step 1
  -p["bPS"] * c(rep(abs(v * 1.5 - sPref)^p["aPS"], 11), 
                rep(abs(v - sPref)^p["aPS"],11),
                rep(abs(v/2 - sPref)^p["aPS"], 11))  
}


#' Walk-beside Utility
#'
#' @param p Numeric vector with subject parameters.
#' @param WB Numeric vector of distances from cells' centers to closest buddy.
#'
#' @return Numeric vector with walk-beside utilities for each cell.
#' @export
#'
wbUtility_r <- function(p, WB) {
  if (is.null(WB)) {
    return(numeric(33))
  }
  
  # Sum of utilities for buddies
  -apply(sapply(1:dim(WB$dists)[1], function(i) {
    p["bWB"] * WB$buddies["angleDisagree", i] * WB$dists[i, ]^p["aWB"]
  }), 1, sum)
}
