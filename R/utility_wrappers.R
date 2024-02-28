# Wrappers for Utility Functions

#' @rdname baUtility_r
#' @param use Whether the R or C++ implementation is called.
baUtility = function(p, BA, use = 'cpp') {
  if (use == 'r' || (exists('predped_env') && predped_env$use == 'r')) {
    return(m4ma::baUtility_r(p, BA))
  } else {
    if (is.null(BA)) return(numeric(33))
    else return(m4ma::baUtility_rcpp(p["aBA"], p["bBA"], BA, as.integer(names(BA))))
  }
}


#' @rdname caUtility_r
#' @param use Whether the R or C++ implementation is called.
caUtility = function(p, use = 'cpp') {
  if (use == 'r' || (exists('predped_env') && predped_env$use == 'r')) {
    return(m4ma::caUtility_r(p))
  } else {
    return(m4ma::caUtility_rcpp(p["aCA"], p["bCA"], p['bCAlr']))
  }
}


#' @rdname flUtility_r
#' @param use Whether the R or C++ implementation is called.
flUtility = function(p, FL, use = 'cpp') {
  if (use == 'r' || (exists('predped_env') && predped_env$use == 'r')) {
    return(m4ma::flUtility_r(p, FL))
  } else {
    if (is.null(FL)) return(numeric(33))
    else return(m4ma::flUtility_rcpp(
        p["aFL"], p["bFL"], p['dFL'], FL$leaders, FL$dists
      )
    )
  }
}


#' @rdname gaUtility_r
#' @param use Whether the R or C++ implementation is called.
gaUtility = function(p, GA, use = 'cpp') {
  if (use == 'r' || (exists('predped_env') && predped_env$use == 'r')) {
    return(m4ma::gaUtility_r(p, GA))
  } else {
    return(m4ma::gaUtility_rcpp(p["bGA"], p["aGA"], GA))
  }
}


#' @rdname idUtility_r
#' @param use Whether the R or C++ implementation is called.
idUtility = function(p, n, ID, ok, group, use = 'cpp') {
  if (use == 'r' || (exists('predped_env') && predped_env$use == 'r')) {
    return(m4ma::idUtility_r(p, n, ID, ok, group))
  } else {
    if (is.null(ID) | !any(ok)) {  
       return(as.vector(ifelse(ok, 0, -Inf)))
    }
    namesInGroup <- names(group[-n][group[-n] == group[n]])
    inGroup <- dimnames(ID)[[1]] %in% namesInGroup
    ok <- ok & apply(ID, 2, function(x) {
      all(x > 0)
    })
    u <- ifelse(ok, 0, -Inf)
    return(m4ma::idUtility_rcpp(p["bID"], p["dID"], p['aID'], inGroup, ok, ID, u))
  }
}


#' @rdname psUtility_r
#' @param use Whether the R or C++ implementation is called.
psUtility = function(p, v, d, use = 'cpp') {
  if (use == 'r' || (exists('predped_env') && predped_env$use == 'r')) {
    return(m4ma::psUtility_r(p, v, d))
  } else {
    return(m4ma::psUtility_rcpp(p["aPS"], p["bPS"], p['sPref'], p['sSlow'], v, d))
  }
}


#' @rdname wbUtility_r
#' @param use Whether the R or C++ implementation is called.
wbUtility = function(p, WB, use = 'cpp') {
  if (use == 'r' || (exists('predped_env') && predped_env$use == 'r')) {
    return(m4ma::wbUtility_r(p, WB))
  } else {
    if (is.null(WB)) return(numeric(33))
    else return(m4ma::wbUtility_rcpp(p["aWB"], p["bWB"], WB$buddies, WB$dists))
  }
}


#' @rdname pCNL_r
#' @param use Whether the R or C++ implementation is called.
pCNL = function(cell, V, muM = rep(1, length(nests)),
                nests, alpha, mu, cellNest, use = 'cpp') {
  if (use == 'r' || (exists('predped_env') && predped_env$use == 'r')) {
    return(m4ma::pCNL_r(cell, V, muM = rep(1, length(nests)),
                        nests, alpha, mu, cellNest))
  } else {
    return(m4ma::pcnl_rcpp(cell_nest[cell, ], V, rep(1, length(nests)),
                           nests, alpha, mu))
  }
}