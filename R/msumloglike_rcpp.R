
get_mum <- function(p) {
  1 / (1 - p[c("Central", "NonCentral", "acc", "const", "dec")])
}

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

like_states_rcpp <- function(trace, p, nests, alpha, cell_nest) {
  
  lhoods <- lapply(trace, function(element) {
    iteration <- sapply(1:length(element$v), function(n) {
      m4ma::like_state_rcpp(element, p, n, nests, alpha, cell_nest)
    })
    return(iteration)
  })
  
  return(lhoods)
}

msumlogLike_rcpp <- function(p, prepped_trace, minLike = 1e-10, 
                             mult = -1) {
  
  imat <- suppressWarnings(matrix(1:length(prepped_trace), ncol = 1))
  # Remove duplicated assignments
  imat[duplicated(as.vector(imat))] <- NA 
  
  # Create lists of prepped_trace, collections of iterations to spread over cores
  prepped_traces <- apply(imat, 2, function(x) {
    prepped_trace[x[!is.na(x)]]
  })
  
  out = lapply(prepped_traces, m4ma::like_states_rcpp, p = p,
               nests = attr(prepped_trace, "nests"),
               alpha = attr(prepped_trace, "alpha"),
               cell_nest = m4ma::get_cell_nest())
  
  return(mult * sum(log(pmax(unlist(out), minLike))))
}