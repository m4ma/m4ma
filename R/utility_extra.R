
#' @rdname destinationAngle_rcpp
destinationAngle_r <- function(a, p1, P1,
                             angles = c(72.5, 50, 32.5, 20, 10, 0, 350, 340, 
                                        327.5, 310, 287.5)) {
  sapply((angles + a) %% 360, minAngle, a2 = angle2(p1, P1))
}


#' @rdname predClose_rcpp
predClose_r <- function(n, p1, a1, p2, r, centres, p_pred, objects) {
  if (dim(p_pred)[1] == 1) {
    return(NULL)
  }
  
  # remove self and pedestrians you cant see 
  occluded <- c(1:dim(p2)[1])[-n][!seesMany(p1, ps = p2[-n, , drop = FALSE], 
                                            objects)]
  p_pred <- p_pred[-c(n,occluded), , drop = FALSE]
  
  if (dim(p_pred)[1] <= 1) {
    return(NULL)
  }
  p2 <- p2[-c(n, occluded), , drop = FALSE]
  
  # Peds in field of vision now and when moved
  inFront <- (minAngle(a1, angle2(p1, p2)) < 85) & 
    (minAngle(a1, angle2(p1, p_pred)) < 85)
  
  if (!any(inFront)) {
    return(NULL)
  }
  
  # Distance to predicted positions
  d <- matrix(apply(centres, 1, dist1, p2 = p_pred[inFront, , drop = FALSE]),
              nrow = sum(inFront), 
              dimnames = list(names(inFront[inFront]), NULL))
  d <- d - (r[n] + r[names(inFront)[inFront]]) # subtract bodies
  d
}


#' @rdname eObjects_rcpp
eObjects_r <- function(p1, p2, r) {
  d <- dist1(p1, p2)
  a12 <- angle2(p1, p2)
  r <- rep(r, length.out = dim(p2)[1])  # in case a single value
  theta <- atan(r / d) * 180 / pi
  ac <- t(as.vector(p1) + t(aTOd((a12 + theta) %% 360) * d))  # anti-clockwise
  cw <- t(as.vector(p1) + t(aTOd((a12 - theta) %% 360) * d))  # clockwise
  # points(acw,pch=16)
  # points(cw,pch=16)
  array(c(ac, cw), dim = c(dim(ac), 2), 
        dimnames = list(names(d), dimnames(ac)[[2]], ends = c("ac", "cw")))
}


#' @rdname iCones_rcpp
iCones_r <- function(p1, a, p2, r, objects) {
  
  # One end in, one out, fill in extreme cone for out
  fix <- function(x) {
    if (all(is.na(x))) {
      return(NA)
    }
    if (all(is.na(x) == c(FALSE, TRUE))) {
      
      x <- c(x[1],11)  # cw out
    }
    if (all(is.na(x) == c(TRUE, FALSE))) {
      x <- c(1,x[2])  # ac out
    }
    if (length(x) > 1) {
      x <- x[1]:x[2]  # fill in intervening
    }
    
    return(x)
  }
  
  if (dim(p2)[1] == 0) {
    return(NULL)
  }
  ends <- eObjects_r(p1, p2, r)
  endCones <- apply(ends, 3, function(x) {
    Iangle_r(p1, a, x)
  })
  if (dim(ends)[1] == 1) { 
    endCones <- matrix(endCones, nrow = 1, 
                       dimnames = list(dimnames(p2)[1], NULL))
  }
  
  cList <- vector(mode = "list", length = dim(endCones)[1])
  names(cList) <- dimnames(endCones)[[1]]
  for (i in 1:length(cList)) {
    cList[[i]] <- fix(endCones[i, ])
  }
  cList <- cList[unlist(lapply(cList, function(x) {
    !any(is.na(x))
  }))]
  if (length(cList) == 0) {
    return(NULL)
  }
  
  # List of candidate cones for each participant in front
  cList <- cList[!unlist(lapply(cList, is.null))]
  if (length(cList) == 0) {
    return(NULL)
  }
  
  # End of unit line from p1 in direction of each cone
  coneLineEnds <- c_vd(1:11, as.vector(p1), rep(1, 11), a)
  print(cList)
  # Distances to objects in cones
  cDist <- vector(mode = "list", length = length(cList)) 
  for (i in names(cList)) {  # remove when cant see
    if (length(cList[[i]]) == 1) {
      if (!seesGoal(p1, p2[i, , drop = FALSE], objects)) {
        cList[[i]] <- numeric(0)
      } else {
        cDist[[i]] <- dist1(p1, p2[i, , drop = FALSE])
      }
    } else {  # more than one cone to check
      for (j in 1:length(cList[[i]])) {  # remove cant see
        # Intersection of cone line and end
        P_n <- matrix(line.line.intersection(p1, coneLineEnds[cList[[i]][j], ], 
                                             ends[i, , 1], ends[i, , 2]), 
                      nrow = 1)
        if (!seesGoal(p1, P_n, objects)) {
          cList[[i]][j] <- NA
          print(TRUE)
        } else {
          if (j == 1) {
            cDist[[i]] <- dist1(p1,P_n) 
          } else {
            cDist[[i]] <- c(cDist[[i]], dist1(p1, P_n))
          }
        }
      }
      cList[[i]] <- cList[[i]][!is.na(cList[[i]])]
    }
  }
  cList <- cList[!unlist(lapply(cList, function(x) {
    length(x) == 0
  }))]
  if (length(cList) == 0) {
    return(NULL)
  }
  cDist <- cDist[!unlist(lapply(cDist, is.null))]
  for (i in 1:length(cList)) {
    names(cDist[[i]]) <- cList[[i]]
  }
  outCones <- out <- numeric(0)
  for (i in 1:11) {
    d <- unlist(lapply(cDist, function(x) {
      x[names(x) == i]
    }))
    if (length(d) > 0) {
      outCones <- c(outCones, i)
      out <- c(out, min(d))
    }
  }
  names(out) <- outCones
  out
}


#' @rdname iCones2Cells_rcpp
iCones2Cells_r <- function(iC, v, vels = c(1.5, 1, .5)) {
  out <- rep(iC, times = 3) - rep(scaleVel(v) * vels, each = length(iC))
  names(out) <- rep(as.numeric(names(iC)), times = 3) + 
    rep(c(0, 11, 22), each = length(iC))
  return(out)
}


#' @rdname blockedAngle_rcpp
blockedAngle_r <- function(n, state, p_pred, objects) {
  iC <- iCones_r(p1 = state$p[n, , drop = FALSE], a = state$a[n], 
               p2 = p_pred[-n, , drop = FALSE], r = state$r, objects)
  iCones2Cells_r(iC, state$v[n])
}


#' @rdname getLeaders_rcpp
getLeaders_r <- function(n, state, centres, objects, onlyGroup = FALSE, 
                       preferGroup = TRUE, pickBest = FALSE) {
  p1 <- state$p[n, , drop = FALSE]
  a1 <- state$a[n] 
  v1 <- state$v[n] 
  
  # Remove peds cant see
  occluded <- c(n, c(1:length(state$v))[-n][
    !seesMany(p1, ps = state$p[-n, , drop = FALSE], objects)])
  a2 <- state$a[-occluded]
  p2 <- state$p[-occluded, , drop = FALSE]
  
  if (dim(p2)[1] == 0) {
    return(NULL)
  }
  if (!is.list(state$P)) {
    P1 <- state$P[n, 1:2, drop = FALSE] 
  } else {
    P1 <- state$P[[n]][attr(state$P[[n]], "i"), 1:2, drop = FALSE]
  }
  
  I_n <- Iangle_r(p1, a1, p2)  # is in cone
  I_n <- I_n[!is.na(I_n)]
  if (length(I_n) == 0) {
    return(NULL)
  }
  
  # Subset in rings 1-3
  ring <- as.numeric(
    as.character(cut(dist1_r(p1, p2[names(I_n), , drop = FALSE]), 
                     c(0, scaleVel(v1) * c(.5, 1, 5)), 
                     c("3", "2", "1"))))
  names(ring) <- names(I_n)
  ring <- ring[!is.na(ring)]
  if (length(ring) == 0) {
    return(NULL)
  }
  
  candidates <- I_n[names(ring)] + 11 * (ring - 1)
  inGroup <- state$group[names(candidates)] == state$group[n]
  if (onlyGroup) {
    if (!any(inGroup)) {
      return(NULL) 
    } else {
      candidates <- candidates[inGroup]
    }
  } else if (preferGroup & any(inGroup)) {
    candidates <- candidates[inGroup]
  }
  
  # Difference in angle between leader heading and destination
  angles <- sapply(a2[names(candidates)], minAngle_r, a2 = Dn_r(p1, P1))
  ok <- angles < 90
  if (!any(ok)) {
    return(NULL)
  }
  candidates <- candidates[ok]
  angles <- angles[ok]
  if (!any(duplicated(candidates))) {
    print(1)
    leaders <- candidates 
  } else {
    print(2)
    leaders <- unique(candidates)
    for (i in 1:length(leaders)) { 
      names(leaders)[i] <- names(candidates)[candidates == leaders[i]][
        which.min(angles[candidates == leaders[i]])] 
    }
    angles <- angles[names(leaders)]
  }
  
  # Distances to leader cells
  d <- array(dim = c(length(leaders), 33), 
             dimnames = list(names(leaders), 1:33))
  for (i in 1:length(leaders)) {
    d[i, ] <- dist1(centres[leaders[i], ], centres)
  }
  if (pickBest) {
    best <- which.min(angles) 
  } else {
    best <- 1:length(angles)
  }
  
  list(dists = d[best, , drop = FALSE], 
       leaders = rbind(cell = leaders, 
                       angleDisagree = (angles[names(leaders)]) / 90, 
                       inGroup = inGroup[names(leaders)])[, best, drop = FALSE])
}


#' @rdname getBuddy_rcpp
getBuddy_r <- function(n, group, a, p_pred, centres, objects, pickBest = FALSE, 
                     state) {
  # Remove peds cant see and self
  occluded <- c(n, c(1:length(state$v))[-n][
    !seesMany(state$p[n, , drop = FALSE], ps = state$p[-n, , drop = FALSE], 
              objects)])
  p_pred <- p_pred[-occluded, , drop = FALSE]
  if (dim(p_pred)[1] == 0) {
    return(NULL)
  }
  
  # NB: Not checking in front as assume you know where your group members are.
  #     If they are behind this formulation will tend to slow you down so they 
  #     can catch up.
  inGroup <- group[-occluded] == group[n]
  p_pred <- p_pred[inGroup, , drop = FALSE]
  nped <- dim(p_pred)[1]
  if (nped == 0) {
    return(NULL)
  }
  
  # Potential buddy matrix of difference between cone angles (rows) 
  # and heading angle for potential buddy (columns)
  headingDifference <- sapply(a[row.names(p_pred)], headingAngle, a1 = a[n])
  
  # Most parallel cones for each potential companion
  parallelCone <- apply(headingDifference, 2, which.min)
  
  # Distances from predicted buddy position to acc, const and dec rings in each 
  # cone
  d <- sapply(1:nped, function(x) {
    dist1(p_pred[names(parallelCone)[x], ],
          centres[c(parallelCone[x], parallelCone[x] + 11, 
                    parallelCone[x] + 22), ])
  })
  
  # Ring closest to buddies
  ring <- apply(d, 2, which.min)
  
  # Cell closest to buddy
  cell <- cbind(parallelCone, parallelCone + 11, 
                parallelCone + 22)[cbind(1:nped, ring)]
  names(cell) <- names(parallelCone)
  
  # Heading difference for each buddy
  angleDisagree <- (headingDifference[cbind(parallelCone, 
                                            1:length(parallelCone))]) / 90
  
  # Distances to buddies cells
  d <- array(dim = c(nped, 33), dimnames = list(names(angleDisagree), 1:33))
  for (i in 1:nped) {
    d[i,] <- dist1(centres[cell[i], ], centres)
  }
  if (pickBest) {
    best <- which.min(angleDisagree) 
  } else {
    best <- 1:nped
  }
  list(buddies = rbind(cell, angleDisagree)[, best, drop = FALSE],
       dists = d[best, , drop = FALSE])
}