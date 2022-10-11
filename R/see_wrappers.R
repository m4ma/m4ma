# Wrappers for See Functions

#' @rdname line_line_intersection_rcpp
#' @param ... Arguments passed to the function implementation.
#' @param use Whether the R or C++ implementation is called.
line.line.intersection = wrapper('line_line_intersection')


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