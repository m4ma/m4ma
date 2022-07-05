baUtility <- function(p, BA) {
  if (is.null(BA)) {
    return(rep(0,33))
  }
  cells <- rep(0, 33)
  cells[as.numeric(names(BA))] <- p["bBA"] / (pmax(BA, 0)^p["aBA"])
  return(-cells)
}