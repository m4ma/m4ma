# Wrappers for See Functions

#' @rdname line_line_intersection_rcpp
#' @param use Whether the R or C++ implementation is called.
line.line.intersection = function(P1, P2, P3, P4, interior.only = FALSE,
                                  use = 'cpp') {
  if (use == 'r' || (exists('predped_env') && predped_env$use == 'r')) {
    return(m4ma::line_line_intersection_r(P1, P2, P3, P4, interior.only))
  } else {
    return(m4ma::line_line_intersection_rcpp(P1, P2, P3, P4,
                                             interior_only = interior.only))
  }
}


#' @rdname seesGoal_rcpp
#' @param ... Arguments passed to the function implementation.
#' @param use Whether the R or C++ implementation is called.
seesGoal = wrapper('seesGoal')


#' @rdname seesCurrentGoal_rcpp
#' @param ... Arguments passed to the function implementation.
#' @param use Whether the R or C++ implementation is called.
seesCurrentGoal = wrapper('seesCurrentGoal')


#' @rdname seesMany_rcpp
#' @param ... Arguments passed to the function implementation.
#' @param use Whether the R or C++ implementation is called.
seesMany = wrapper('seesMany')


#' @rdname seesGoalOK_rcpp
#' @param ... Arguments passed to the function implementation.
#' @param use Whether the R or C++ implementation is called.
seesGoalOK = wrapper('seesGoalOK')