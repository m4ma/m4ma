
#' Nest Association
#'
#' @param p Numeric vector with subject parameters.
#'
#' @return Numeric vector with nest association parameters.
#' @export
#'
get_mum <- function(p) {
  1 / (1 - p[c("Central", "NonCentral", "acc", "const", "dec")])
}

#' Cell Nest Association
#'
#' @return Integer matrix indicating which cells belong to which nests. Indices start at zero to match C++ indexing.
#' @export
#'
get_cell_nest = function() {
  cell_nest_rcpp = cbind(
    t(matrix(c(c(1, 5),           # cell 0: Central, dec
               rep(2:3, 5),       # cell 1-5: NonCentral, acc
               c(1, 3),           # cell 6: Central, acc
               rep(2:3, 5),       # cell 7-11: NonCentral, acc
               rep(c(2, 4), 5),   # cell 12-16: NonCentral, const
               c(1, 4),           # cell 17: Central, const
               rep(c(2, 4), 5),   # cell 18-22: NonCentral, const
               rep(c(2, 5), 5),   # cell 23-27: NonCentral, dec
               c(1, 5),           # cell 28: Central, dec
               rep(c(2, 5), 5)),  # cell 29-33: NonCentral, dec
             nrow = 2)),
    # Index cell in NonCentral/Central and acc/const/dec nest
    cbind(c(1, 1:5, 2, 6:15, 3, 16:25, 4, 26:30), 
          c(1, 1:11, 1:11, 2:12))) - 1 # Index needs to start from zero for Rcpp
  
  return(cell_nest_rcpp)
}

#' State Likelihood
#'
#' @param state List of state variables.
#' @param p Numeric vector of subject parameters.
#' @param n Integer scalar indexing the subject in the state.
#' @param nests List of vectors with utility indices.
#' @param alpha List of vectors with alpha values.
#' @param cell_nest Integer matrix associating the cells with the nests.
#'
#' @return Probability of the state given the parameters.
#' @export
#'
like_state_rcpp <- function(state, p, n, nests, alpha, cell_nest) {
  p_n <- p[names(state$v)[n], ]
  
  u <- m4ma::utility_rcpp(p_n, n, state$v[n], state$d[n], state$GA[[n]],
                          state$BA[[n]], state$ID[[n]], state$FL[[n]],
                          state$WB[[n]], state$ok[[n]], state$group)
  cell <- state$cell[n]
  
  prob <- m4ma::pcnl_rcpp(cell_nest[cell + 1, ], u, mum = m4ma::get_mum(p_n),
                          nests, alpha, mu = 1)
  
  return(prob)
}

#' Trace Likelihood
#'
#' @param trace List of lists with state variables.
#' @param p Numeric matrix with subject parameters.
#' @param nests List of vectors with utility indices.
#' @param alpha List of vectors with alpha values.
#' @param cell_nest Integer matrix associating the cells with the nests.
#'
#' @return List with a numeric vector for each state containing the state likelihoods.
#' @export
#'
like_states_rcpp <- function(trace, p, nests, alpha, cell_nest) {
  
  lhoods <- lapply(trace, function(element) {
    iteration <- sapply(1:length(element$v), function(n) {
      m4ma::like_state_rcpp(element, p, n, nests, alpha, cell_nest)
    })
    return(iteration)
  })
  
  return(lhoods)
}

#' Log-likelihood Sum
#'
#' @param trace List of lists with state variables.
#' @param p Numeric matrix with subject parameters.
#' @param minLike Numeric scalar indicating minimum likelihood value.
#' @param mult Numeric scalar indicating likelihood sum multiplicator.
#'
#' @return Numeric scalar indicating log-likelihood sum of trace given parameters.
#' @export
#'
msumlogLike_rcpp <- function(trace, p, minLike = 1e-10, mult = -1) {
  
  out = m4ma::like_states_rcpp(
    p = p,
    trace = trace,
    nests = attr(trace, "nests"),
    alpha = attr(trace, "alpha"),
    cell_nest = m4ma::get_cell_nest()
  )
  
  return(mult * sum(log(pmax(unlist(out), minLike))))
}