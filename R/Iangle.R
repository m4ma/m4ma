#' Iangle
#' 
#' Which angle cone (1..11, NA means outside of view) is p2 in relative to p1 heading at angel a1 
#' @param p1 Numeric matrix
#' @param p2 Numeric matrix
#' @param a1 Numeric vector - angles in degrees between 0 and 360
#' @param border Numeric vector - perimeter boundaries
#' @return Numeric Vector of indices whereby the angle is within border's bins
Iangle <- function(p1, a1, p2, border = c(-85, -60, -40, -25, -15, -5, 5, 15, 
                                          25, 40, 60, 85)) {
  
  tomp <- function(x) { # -180 ... 180
    neg <- x < -180
    x[neg] <- 360 + x[neg]
    
    pos <- x > 180
    x[pos] <- x[pos] - 360
    
    x
  }
  
  a <- (angle2(p1, p2) - a1)
  
  out <- 12 - .bincode(tomp(a), border)
  names(out) <- row.names(p2)
  out
} 