#' Total Utility of Alternatives
#'
#' @param p Vector with subject parameters.
#' @param v Scalar indicating the current speed..
#' @param d Scalar indicating the distance to next goal.
#' @param ga Vector with angles to next goal.
#'
#' @return Vector with total utility for each alternative.
#' @export
#'
utility_rcpp = function(p, n, v, d, ga, ba, id, fl, wb, ok, group) {
  ps_utility = m4ma::psUtility_rcpp(p['aPS'], p['sPref'], p['sSlow'], p['bPS'], v, d)
  ga_utility = m4ma::gaUtility_rcpp(p['bGA'], p['aGA'], ga)
  ca_utility = m4ma::caUtility_rcpp(p['aCA'], p['bCA'], p['bCAlr'])
  
  if (!is.null(ba))
    ba_utility = m4ma::baUtility_rcpp(p['aBA'], p['bBA'], ba, as.integer(names(ba)) - 1)
  else
    ba_utility = numeric(33)
  
  id_utility = m4ma::idUtility(p, n, id, ok, group)
  
  if (!is.null(fl))
    fl_utility = m4ma::flUtility_rcpp(p['aFL'], p['bFL'], p['dFL'], fl$leaders, fl$dists)
  else
    fl_utility = numeric(33)
  
  wb_utility = m4ma::wbUtility(p, wb)
  
  total_utility = ps_utility + ga_utility + ca_utility + ba_utility + 
    id_utility + fl_utility + wb_utility
  
  return(c(-p['bS'], total_utility) / p['rU'])
}
