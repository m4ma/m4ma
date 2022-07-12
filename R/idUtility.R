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
idUtility <- function(p, n, ID, ok, group) {
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