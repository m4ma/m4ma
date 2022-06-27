caUtility <- function(p, angles = c(10, 20, 32.5, 50, 72.5) / 90) {
  ap <- angles^p["aCA"]
  -rep(c(rep(p["bCA"] * p["bCAlr"], 5), 1, rep(p["bCA"] / p["bCAlr"], 5)) * 
         c(ap[5:1], 0, ap), times = 3)
}
