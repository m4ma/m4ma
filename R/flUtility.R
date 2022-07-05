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