# Anti-clockwise angle from p1 as origin to p2 (x,y pairs matrices) 0 to <360
angle2 <- function(p1, p2) {
  #round((atan2(p2[, 2] - p1[, 2], p2[, 1] - p1[, 1]) * 180 / pi), 10) %% 360
  (atan2(p2[, 2] - p1[, 2], p2[, 1] - p1[, 1]) * 180 / pi) %% 360
}
