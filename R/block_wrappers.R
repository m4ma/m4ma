# Wrappers for block functions

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


#' @rdname object2lines_rcpp
#' @param ... Arguments passed to the function implementation.
#' @param use Whether the R or C++ implementation is called.
object2lines = wrapper('object2lines')


#' @rdname bodyObjectOverlap_rcpp
#' @param ... Arguments passed to the function implementation.
#' @param use Whether the R or C++ implementation is called.
bodyObjectOverlap = wrapper('bodyObjectOverlap')


#' @rdname bodyObjectOK_rcpp
#' @param ... Arguments passed to the function implementation.
#' @param use Whether the R or C++ implementation is called.
bodyObjectOK = wrapper('bodyObjectOK')