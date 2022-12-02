# Wrappers for Utility Helper Functions

#' @rdname destinationAngle_rcpp
#' @param ... Arguments passed to the function implementation.
#' @param use Whether the R or C++ implementation is called.
destinationAngle = wrapper('destinationAngle')


#' @rdname predClose_rcpp
#' @param ... Arguments passed to the function implementation.
#' @param use Whether the R or C++ implementation is called.
predClose = wrapper('predClose')


#' @rdname eObjects_rcpp
#' @param ... Arguments passed to the function implementation.
#' @param use Whether the R or C++ implementation is called.
eObjects = wrapper('eObjects')


#' @rdname iCones_rcpp
#' @param ... Arguments passed to the function implementation.
#' @param use Whether the R or C++ implementation is called.
iCones = wrapper('iCones')


#' @rdname iCones2Cells_rcpp
#' @param ... Arguments passed to the function implementation.
#' @param use Whether the R or C++ implementation is called.
iCones2Cells = wrapper('iCones2Cells')


#' @rdname blockedAngle_rcpp
#' @param ... Arguments passed to the function implementation.
#' @param use Whether the R or C++ implementation is called.
blockedAngle = wrapper('blockedAngle')


#' @rdname getLeaders_rcpp
#' @param ... Arguments passed to the function implementation.
#' @param use Whether the R or C++ implementation is called.
getLeaders = wrapper('getLeaders')


#' @rdname getBuddy_rcpp
#' @param use Whether the R or C++ implementation is called.
getBuddy = function(n, group, a, p_pred, centres, objects, pickBest = FALSE, 
                    state, use = 'cpp') {
  if (use == 'r' || (exists('predped_env') && predped_env$use == 'r')) {
    return(m4ma::getBuddy_r(n, group, a, p_pred, centres, objects,
                            pickBest, state))
  } else {
    return(m4ma::getBuddy_rcpp(n, group, a, p_pred, centres, objects,
                               pickBest, state))
  }
}