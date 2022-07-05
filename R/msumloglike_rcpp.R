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