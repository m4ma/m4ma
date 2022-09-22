# Wrappers for Geometry Functions

#' Execute R or C++ Implementation
#' 
#' This is a wrapper function that takes a function name and returns a function
#' that calls the function with an additional argument `use` to switch between 
#' an R and C++ version.
#'
#' @param name The name of the function to be wrapped.
#'
#' @return The wrapped function with the additional `use` argument to switch
#' implementations.
#' 
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
#' @param ... Arguments passed to the function implementation.
#' @param use Whether the R or C++ implementation is called.
dist = wrapper('dist')


#' @rdname dist1_rcpp
#' @param ... Arguments passed to the function implementation.
#' @param use Whether the R or C++ implementation is called.
dist1 = wrapper('dist1')


#' @rdname angle2s_rcpp
#' @param ... Arguments passed to the function implementation.
#' @param use Whether the R or C++ implementation is called.
angle2s = wrapper('angle2s')


#' @rdname angle2s_rcpp
#' @param ... Arguments passed to the function implementation.
#' @param use Whether the R or C++ implementation is called.
angle2 = wrapper('angle2')


#' @rdname aTOd_rcpp
#' @param ... Arguments passed to the function implementation.
#' @param use Whether the R or C++ implementation is called.
aTOd = wrapper('aTOd')


#' @rdname Iangle_rcpp
#' @param ... Arguments passed to the function implementation.
#' @param use Whether the R or C++ implementation is called.
Iangle = wrapper('Iangle')


#' @rdname angle2s_rcpp
#' @param ... Arguments passed to the function implementation.
#' @param use Whether the R or C++ implementation is called.
Dn = wrapper('Dn')


#' @rdname minAngle_rcpp
#' @param ... Arguments passed to the function implementation.
#' @param use Whether the R or C++ implementation is called.
minAngle = wrapper('minAngle')


#' @rdname headingAngle_rcpp
#' @param ... Arguments passed to the function implementation.
#' @param use Whether the R or C++ implementation is called.
headingAngle = wrapper('headingAngle')


#' @rdname scaleVel_rcpp
#' @param ... Arguments passed to the function implementation.
#' @param use Whether the R or C++ implementation is called.
scaleVel = wrapper('scaleVel')


#' @rdname c_vd_rcpp
#' @param ... Arguments passed to the function implementation.
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
#' @param ... Arguments passed to the function implementation.
#' @param use Whether the R or C++ implementation is called.
coneNum = wrapper('coneNum')


#' @rdname ringNum_rcpp
#' @param ... Arguments passed to the function implementation.
#' @param use Whether the R or C++ implementation is called.
ringNum = wrapper('ringNum')