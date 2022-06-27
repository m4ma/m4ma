# Goal angle utility, GA = goal angle (degrees/90) for 11 cones  
gaUtility <- function(p, GA) {
  -rep(p["bGA"] * GA^p["aGA"], times = 3)
}

