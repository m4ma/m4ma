# Angle to destination for all pedestrians
Dn <- function(p_n, P_n) {
  out <- numeric(dim(p_n)[1])
  for (i in 1:dim(p_n)[1]) {
    out[i] <- angle2(p_n[i, , drop = FALSE], P_n[i, , drop = FALSE])
  }
  out
}