# Shortest angle anti-clockwise from p1 as orign to p2 (> -180 to 180)
angle2s <- function(p1, p2) {
  round((atan2(p2[, 2] - p1[, 2], p2[, 1] - p1[, 1]) * 180 / pi), 10)
}
