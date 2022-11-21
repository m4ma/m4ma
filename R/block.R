
#' @rdname object2lines_rcpp
object2lines_r <- function(o) {
  array(c(o$x[1], o$y[1], o$x[1], o$y[2], o$x[1], o$y[1], o$x[2], o$y[1],
          o$x[2], o$y[1], o$x[2], o$y[2], o$x[1], o$y[2], o$x[2], o$y[2]),  
        dim = c(2, 4, 2), dimnames = list(c("x", "y"), 
                                          c("L1", "L2", "L3", "L4"), 
                                          c("P1", "P2")))
}


#' @rdname bodyObjectOverlap_rcpp
bodyObjectOverlap_r <- function(oL, r, okCentres) {
  # Right angles to each object line
  a <- unique((angle2_r(p1 = t(oL[, , "P1"]), p2 = t(oL[, , "P2"])) + 90) %% 180)
  # dx and dy to move along line
  x <- r * sin(a * pi / 180)
  y <- r * cos(a * pi / 180)
  
  # For each centre check for overlap
  apply(okCentres, 1, function(p) { 
    # Lines 
    segments <- array(c(p - rbind(x, y), p + rbind(x, y)), 
                      dim = c(2, length(x), 2),
                      dimnames = list(c("x", "y"), NULL, c("P1", "P2")))
    out <- FALSE 
    for (j in 1:dim(segments)[2]) {
      out <- out | any(sapply(1:(dim(oL)[2]), function(i){
        all(is.finite(line_line_intersection_r(segments[,j,"P1"], 
                                               segments[,j,"P2"], 
                                               oL[,i,"P1"], oL[,i,"P2"], TRUE)
        )) 
      }))
    }
    out
  })
}


#' @rdname bodyObjectOK_rcpp
bodyObjectOK_r <- function(r, centres, objects, ok) {
  if (!any(ok)) {
    return(NULL)
  }
  oLines <- lapply(objects, object2lines_r)
  
  # Does it overlap
  out <- !logical(33)
  out[ok] <- apply(matrix(unlist(
    lapply(oLines, bodyObjectOverlap_r, r = r,
           okCentres = centres[as.vector(ok), , drop = FALSE])),
    nrow = sum(ok)), 1, function(x) {
      any(x)
    })
  # If it doesn't overlap it is OK
  matrix(!out, nrow = 11)
}