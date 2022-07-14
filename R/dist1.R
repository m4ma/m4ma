dist1 <- function(p1, p2) {
  sqrt(apply((t(p2) - as.vector(p1))^2, 2, sum))
}