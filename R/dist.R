dist <- function(p1, p2) {
  sqrt(apply((p1 - p2)^2, 1, sum))
}