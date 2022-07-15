# Helper functions for testing

create_rcpp_states = function(states) {
  rcpp_states = lapply(1:length(states$v), function(i) {
    rcpp_state = list(
      v = states$v[i],
      d = states$d[i],
      BA = states$BA[[i]],
      GA = states$GA[[i]],
      ID = states$ID[[i]],
      FL = states$FL[[i]],
      WB = states$WB[[i]],
      ok = states$ok[[i]],
      group = states$group,
      cell = states$cell[i]
    )
    
    return(rcpp_state)
  })
  
  return(rcpp_states)
}

create_rcpp_trace = function(trace) {
  return(lapply(trace, create_rcpp_states))
}
