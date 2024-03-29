
#' Rcpp Trace
#' 
#' Helper function that transforms a trace list simulated by predped into the
#' format required to do likelihood estimation.
#' 
#' @param trace List of states in the simulation.
#' @param elements Character string indicating whether to return a trace with 
#' iterations or subjects as elements. Can be either \code{"iterations"} or 
#' \code{"subjects"}.
#' @return List of simulation iterations containing subject states or subjects
#' containing iteration states.
create_rcpp_trace <- function (trace, elements = "iterations") 
{
  pMat <- attr(trace,"pMat")
    if (elements == "iterations") {
      rcpp_trace = lapply(trace, function(states) {
            pn <- which(dimnames(pMat)[[1]] %in% names(states$v))  
            rcpp_states = lapply(1:length(states$v), function(i) {
                if (is.null(states$WB[[i]][[1]])) 
                  WB <- NULL
                else WB <- states$WB[[i]]
                if (length(states$BA[[i]]) == 0) 
                  BA <- NULL
                else BA <- pmax(states$BA[[i]], 0)
                ok <- states$ok[[i]]
                if (!is.null(states$ID[[i]])) {
                  group <- states$group
                  namesInGroup <- names(group[-i][group[-i] == group[i]])
                  inGroup <- dimnames(states$ID[[i]])[[1]] %in% namesInGroup
                  ok <- ok & apply(states$ID[[i]], 2, function(x) {
                    all(x > 0)
                  })
                } else {
                  inGroup <- FALSE
                }
                rcpp_state = list(v = states$v[i], d = states$d[i], 
                  BA = BA, GA = states$GA[[i]], ID = states$ID[[i]],
                  inGroup = inGroup, IDPreUtility = as.vector(ifelse(ok, 0, -Inf)),
                  okID = ok, FL = states$FL[[i]], WB = WB, ok = states$ok[[i]], 
                  group = states$group, cell = states$cell[i], pn = pn)
                return(rcpp_state)
            })
            return(rcpp_states);
        })
        attributes(rcpp_trace) = attributes(trace)
    } else {
        subject_names = attr(trace, "subject_names")
        rcpp_trace = vector(mode = "list", length = length(subject_names))
        names(rcpp_trace) = subject_names
        for (s_name in subject_names) {
            tmp = lapply(trace, function(state, s) {
                n = which(row.names(state$p) == s)
                if (length(n) > 0) {
                  if (is.null(state$WB[[n]][[1]])) 
                    WB <- NULL
                  else WB <- state$WB[[n]]
                  if (length(state$BA[[n]]) == 0) 
                    BA <- NULL
                  else BA <- pmax(state$BA[[n]], 0)
                  ok <- state$ok[[n]]
                  if (!is.null(state$ID[[n]])) {
                    group <- state$group
                    namesInGroup <- names(group[-n][group[-n] == group[n]])
                    inGroup <- dimnames(state$ID[[n]])[[1]] %in% namesInGroup
                    ok <- ok & apply(state$ID[[n]], 2, function(x) {
                      all(x > 0)
                    })
                  } else {
                    inGroup <- FALSE
                  }
                  rcpp_state <- list(v = state$v[n], d = state$d[n], 
                    group = state$group, cell = state$cell[n], 
                    ok = state$ok[[n]], GA = state$GA[[n]], ID = state$ID[[n]],
                    inGroup = inGroup, IDPreUtility = as.vector(ifelse(ok, 0, -Inf)),
                    okID = ok, BA = BA, FL = state$FL[[n]], n = n, WB = WB)
                }
                else {
                  rcpp_state = NULL
                }
                return(rcpp_state)
            }, s = s_name)
            rcpp_trace[[s_name]] = tmp[!unlist(lapply(tmp, is.null))]
        }
        attributes(rcpp_trace) = attributes(trace)
        attr(rcpp_trace, "elements") = "subjects"
        names(rcpp_trace) <- subject_names
    }
    return(rcpp_trace)
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