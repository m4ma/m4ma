#' @rdname line_line_intersection_rcpp
line_line_intersection_r <- function(P1, P2, P3, P4, interior.only=FALSE) 
{
  P1 <- as.vector(P1)
  P2 <- as.vector(P2)
  P3 <- as.vector(P3)
  P4 <- as.vector(P4)
  
  dx1 <- P1[1] - P2[1]
  dx2 <- P3[1] - P4[1]
  dy1 <- P1[2] - P2[2]
  dy2 <- P3[2] - P4[2]
  
  D <- det(rbind(c(dx1, dy1),
                 c(dx2, dy2)))
  if (D==0) {
    return(c(Inf, Inf))
  }
  D1 <- det(rbind(P1, P2))
  D2 <- det(rbind(P3, P4))
  
  X <- det(rbind(c(D1, dx1),
                 c(D2, dx2)))/D
  Y <- det(rbind(c(D1, dy1),
                 c(D2, dy2)))/D
  
  if (interior.only) {
    ## Compute the fractions of L1 and L2 at which the intersection
    ## occurs
    lambda1 <- -((X-P1[1])*dx1 + (Y-P1[2])*dy1)/(dx1^2 + dy1^2)
    lambda2 <- -((X-P3[1])*dx2 + (Y-P3[2])*dy2)/(dx2^2 + dy2^2)
    if (!((lambda1>0) & (lambda1<1) &
          (lambda2>0) & (lambda2<1))) {
      return(c(NA, NA))
    }
  }
  return(c(X, Y))
}

# Can p see goal P (in any direction), or more generally any point P_n, 
# or is it occluded by objects?
#' @rdname seesGoal_rcpp
seesGoal_r <- function(p_n, P_n, objects) {
  for (i in 1:length(objects)) {
    intersects <- c(
      line.line.intersection(p_n, P_n, c(objects[[i]]$x[1], objects[[i]]$y[1]), 
                             c(objects[[i]]$x[2], objects[[i]]$y[1]), 
                             interior.only = TRUE),
      line.line.intersection(p_n, P_n, c(objects[[i]]$x[1], objects[[i]]$y[1]), 
                             c(objects[[i]]$x[1], objects[[i]]$y[2]), 
                             interior.only = TRUE),
      line.line.intersection(p_n, P_n, c(objects[[i]]$x[1], objects[[i]]$y[2]), 
                             c(objects[[i]]$x[2], objects[[i]]$y[2]), 
                             interior.only = TRUE),
      line.line.intersection(p_n, P_n, c(objects[[i]]$x[2], objects[[i]]$y[1]), 
                             c(objects[[i]]$x[2], objects[[i]]$y[2]), 
                             interior.only = TRUE))
    
    if (any(is.finite(intersects))) {
      return(FALSE)
    }
  }
  return(TRUE)  # if none of the objects intersect the view
}

# seesGoal for current goal
#' @rdname seesCurrentGoal_rcpp
seesCurrentGoal_r <- function(n, state, objects, offset = 0) {
  print(state$P[[n]][attr(state$P[[n]], "i") + offset, 1:2])
  seesGoal(state$p[n, ], state$P[[n]][attr(state$P[[n]], "i") + offset, 1:2],
           objects)
}

# Can position p see points in ps matrix (x,y columns)
#' @rdname seesMany_rcpp
seesMany_r <- function(p1, ps, objects) {
  apply(ps, 1, seesGoal, P_n = p1, objects = objects)  
}

# Boolean indicating if goal is visible from ok cells
#' @rdname seesGoalOK_rcpp
seesGoalOK_r <- function(n, objects, state, centres, ok) {
  if (any(ok)) {
    for (i in c(1:33)[ok]) {
      ok[i] <- seesGoal(centres[i, ], 
                        state$P[[n]][attr(state$P[[n]], "i"), 1:2], 
                        objects) 
    }
  }
  return(ok)
}