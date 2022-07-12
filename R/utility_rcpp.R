#' Total Utility of Cells
#'
#' @param p Numeric vector of subject parameters.
#' @param n Integer scalar indexing the subject in the state.
#' @param v Numeric scalar indicating the current speed.
#' @param d Numeric scalar indicating the distance to next goal.
#' @param ga Numeric vector of angles to next goal.
#' @param ba Numeric vector of distances from each cell to closest pedestrian.
#' @param id Numeric matrix of predicted distances from the subject to other pedestrians in the front.
#' @param fl List of numeric matrices:
#' \describe{
#'   \item{leaders}{Matrix with columns per leader and rows of their normalized angle disagreement and in-group status.}
#'   \item{dists}{Matrix with rows per leader and columns per cell with distances from each cell to chosen cell.}
#' }
#' @param wb Numeric vector of distances from cells' centers to closest buddy.
#' @param ok Logical matrix indicating if cells are blocked.
#' @param group Integer vector with group indices for each pedestrian.
#'
#' @return Numeric vector with total utility for each cell.
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
