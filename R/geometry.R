
#' @rdname dist_rcpp
dist_r <- function(p1, p2) {
  sqrt(apply((p1 - p2)^2, 1, sum))
}


#' @rdname dist1_rcpp
dist1_r <- function(p1, p2) {
  sqrt(apply((t(p2) - as.vector(p1))^2, 2, sum))
}


#' @rdname angle2s_rcpp
angle2s_r <- function(p1, p2) {
  round((atan2(p2[, 2] - p1[, 2], p2[, 1] - p1[, 1]) * 180 / pi), 10)
}


#' @rdname angle2s_rcpp
angle2_r <- function(p1, p2) {
  #round((atan2(p2[, 2] - p1[, 2], p2[, 1] - p1[, 1]) * 180 / pi), 10) %% 360
  (atan2(p2[, 2] - p1[, 2], p2[, 1] - p1[, 1]) * 180 / pi) %% 360
}


#' @rdname aTOd_rcpp
aTOd_r <- function(a) {
  cbind(x = cos(a * pi / 180), y = sin(a * pi / 180))
}


#' @rdname Iangle_rcpp
Iangle_r <- function(p1, a1, p2, border = c(-85, -60, -40, -25, -15, -5, 5, 15, 
                                          25, 40, 60, 85)) {
  
  tomp <- function(x) { # -180 ... 180
    neg <- x < -180
    x[neg] <- 360 + x[neg]
    
    pos <- x > 180
    x[pos] <- x[pos] - 360
    
    x
  }
  
  a <- (angle2_r(p1, p2) - a1)
  
  out <- 12 - .bincode(tomp(a), border)
  names(out) <- row.names(p2)
  out
}


#' @rdname angle2s_rcpp
Dn_r <- function(p_n, P_n) {
  out <- numeric(dim(p_n)[1])
  for (i in 1:dim(p_n)[1]) {
    out[i] <- angle2_r(p_n[i, , drop = FALSE], P_n[i, , drop = FALSE])
  }
  out
}


#' @rdname minAngle_rcpp
minAngle_r <- function(a1,a2) {
  pmin(abs(a1 - a2), abs(pmin(a1, a2) + (360 - pmax(a1, a2))))
}


#' @rdname headingAngle_rcpp
headingAngle_r <- function(a2, a1, angles = c(72.5, 50, 32.5, 20, 10, 0, 350, 
                                            340, 327.5, 310, 287.5)) {
  sapply((angles + a1) %% 360, minAngle_r, a2 = a2)  
}


#' @rdname scaleVel_rcpp
scaleVel_r <- function(v, tStep = 0.5) {
  v * tStep  
}


#' @rdname c_vd_rcpp
c_vd_r <- function(cells, p1, v1, a1,
                 vels = matrix(rep(c(1.5, 1, .5), each = 11), ncol = 3),
                 angles = matrix(rep(c(72.5, 50, 32.5, 20, 10, 0, 350, 340, 
                                       327.5, 310, 287.5), times = 3), 
                                 ncol = 3)) {
  t(p1 + t(scaleVel_r(v1) * vels[cells] * aTOd_r((angles[cells] + a1) %% 360)))
}


#' @rdname coneNum_rcpp
coneNum_r <- function(k) {
  1 + ((k-1) %% 11)
}


#' @rdname ringNum_rcpp
ringNum_r <- function(k) {
  1 + (k-1) %/% 11
}