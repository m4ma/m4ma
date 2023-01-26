
#' @rdname object2lines_rcpp
object2lines_r <- function(o) {
  array(c(o$x[1], o$y[1], o$x[1], o$y[2], o$x[1], o$y[1], o$x[2], o$y[1],
          o$x[2], o$y[1], o$x[2], o$y[2], o$x[1], o$y[2], o$x[2], o$y[2]),  
        dim = c(2, 4, 2), dimnames = list(c("x", "y"), 
                                          c("L1", "L2", "L3", "L4"), 
                                          c("P1", "P2")))
}


#' @rdname bodyObjectOverlap_rcpp
bodyObjectOverlap_r <- function(oL,r, okCentres)
  # No body-object overlap (i.e., no overlap with any line making up the object)
  # boolean vector for each okCentres
{
  nL <- dim(oL)[2]
  nC <- dim(okCentres)[1]
  d21 <- apply(oL,1:2,diff)
  
  d01 <- xy <- dxy <- array(dim=c(2,nL,nC),dimnames=list(c("x","y"),NULL,NULL))
  param <- array(dim=c(nL,nC))
  for (i in 1:nL)  {
    d01[,i,] <- as.vector(t(okCentres)) - as.vector(oL[,i,1])
    param[i,] <- d21[,i]%*%d01[,i,]
  }
  len_sq <- apply(d21^2,2,sum)
  len0 <- len_sq == 0
  param[!len0,] <- param[!len0,]/len_sq[!len0]
  for (j in 1:nC) {
    for (i in 1:nL) {
      if (param[i,j] < 0) xy[,i,j] <- oL[,i,1] else
        if (param[i,j] > 1) xy[,i,j] <- oL[,i,2] else
          xy[,i,j] <- oL[,i,1] + param[i,j]*d21[,i]
      dxy[,i,] <- t(okCentres) - xy[,i,]
    }
  }
  len <- sqrt(apply(dxy^2,2:3,sum))
  apply(r>len,2,any)
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