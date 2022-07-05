wbUtility <- function(p, WB) {
  if (is.null(WB)) {
    return(numeric(33))
  }
  
  # Sum of utilities for buddies
  -apply(sapply(1:dim(WB$dists)[1], function(i) {
    p["bWB"] * WB$buddies["angleDisagree", i] * WB$dists[i, ]^p["aWB"]
  }), 1, sum)
}