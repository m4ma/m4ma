#' Rcpp States
#' 
#' Helper function that transforms a state list simulated by predped into the
#' format required to do likelihood estimation.
#' 
#' @param states List of variables representing a
#' state in the simulation.
#' @return List of subjects states.
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


#' Rcpp Trace
#' 
#' Helper function that transforms a trace list simulated by predped into the
#' format required to do likelihood estimation.
#' 
#' @param trace List of states in the simulation.
#' @return List of simulation iterations containing subject states.
create_rcpp_trace = function(trace) {
  return(lapply(trace, create_rcpp_states))
}


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