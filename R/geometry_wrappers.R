
# This is a wrapper function that takes a function name and returns a function
# that calls the function with an additional argument `use` to switch between
# an R and C++ version.
wrapper = function(name) {
  fun = function(..., use = 'cpp') {
    if (use == 'r' || (exists('predped_env') && predped_env$use == 'r')) {
      return(get(paste0(name, '_r'))(...))
    } else {
      return(get(paste0(name, '_rcpp'))(...))
    }
  }
  
  return(fun)
}


#' @rdname dist_rcpp
dist = wrapper('dist')


#' @rdname dist1_rcpp
dist1 = wrapper('dist1')


#' @rdname angle2s_rcpp
angle2s = wrapper('angle2s')


#' @rdname angle2s_rcpp
angle2 = wrapper('angle2')


#' @rdname aTOd_rcpp
aTOd = wrapper('aTOd')


#' @rdname Iangle_rcpp
Iangle = wrapper('Iangle')


#' @rdname angle2s_rcpp
Dn = wrapper('Dn')


#' @rdname minAngle_rcpp
minAngle = wrapper('minAngle')


#' @rdname headingAngle_rcpp
headingAngle = wrapper('headingAngle')


#' @rdname scaleVel_rcpp
scaleVel = wrapper('scaleVel')


#' @rdname c_vd_rcpp
c_vd = function(...) {
  # This function needs an extra wrap to supply defaults for 'vels' and 'angles' for rcpp function
  fun = wrapper('c_vd')(..., 
                        vels = matrix(rep(c(1.5, 1, 0.5), each = 11), ncol = 3),
                        angles = matrix(
                          rep(
                            c(72.5, 50, 32.5, 20, 10, 0, 350,
                              340, 327.5, 310, 287.5), 
                            times = 3),
                          ncol = 3)
                        )
  return(fun)
}


#' @rdname coneNum_rcpp
coneNum = wrapper('coneNum')


#' @rdname ringNum_rcpp
ringNum = wrapper('ringNum')