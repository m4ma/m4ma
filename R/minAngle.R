# Shortest absolute angle between a1 and a2
minAngle <- function(a1,a2) {
  pmin(abs(a1 - a2), abs(pmin(a1, a2) + (360 - pmax(a1, a2))))
}