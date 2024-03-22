// returns a numeric vector instead of 1xN and Nx1 matrix
#define RCPP_ARMADILLO_RETURN_ROWVEC_AS_VECTOR
#define RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR
#include <RcppArmadillo.h>
#include "utils.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppArmadillo)]]


//' Blocked-angle Utility
//'
//' @param aBA Numeric scalar power parameter.
//' @param bBA Numeric scalar scale parameter.
//' @param BA Numeric vector of distances from cells to closest pedestrians. Must
//' be positive.
//' @param idx_BA Integer vector of cell indices.
//' @param n_cells Integer number of cells for which the utility is calculated. 
//'
//' @return Numeric vector of utilities for each cell.
//' @export
// [[Rcpp::export]]
NumericVector baUtility_rcpp(double aBA, double bBA, NumericVector BA, IntegerVector idx_BA, int n_cells = 33) {
  NumericVector utility (n_cells);
  
  int k = BA.length();
  
  for(int i = 0; i < k; i++) {
    int idx = idx_BA[i];
    utility[idx] = -bBA / std::pow(BA[i], aBA);
  }
  
  return(utility);
}


//' Current-angle Utility
//' 
//' @param aCA Numeric scalar power parameter.
//' @param bCA Numeric scalar scale parameter.
//' @param bCAlr Numeric scalar scale parameter.
//' @param angles Numeric vector with angles for which to calculate utility.
//' @param n_rings Integer number of cell rings for which to calculate utlity.
//' 
//' @return Numeric vector of utilities for each cell.
//' @export
// [[Rcpp::export]]
NumericVector caUtility_rcpp(
    double aCA,
    double bCA,
    double bCAlr,
    NumericVector angles = NumericVector::create(10, 20, 32.5, 50, 72.5) / 90.0,
    int n_rings = 3
) {
  
  // allocate variables
  NumericVector output(2 * angles.length() + 1);
  
  // compute power of angles
  NumericVector ap = pow(angles, aCA);
  
  // fill the output
  int ap_len = ap.length();
  
  NumericVector bpb = rep(bCA * bCAlr, ap_len);
  NumericVector bdb = rep(bCA / bCAlr, ap_len);
  
  for(int i = 0;  i < ap_len; ++i) {
    output[i] = bpb[i] * ap[ap_len - i - 1];
    output[i+ap_len+1] = bdb[i] * ap[i];
  }
  output[ap_len] = 0;
  
  output = -rep(output, n_rings);
  return output;
  
}


//' Follow-leader Utility
//'
//' @param aFL Numeric scalar power parameter.
//' @param bFL Numeric scalar scale parameter.
//' @param dFL Numeric scalar direction parameter.
//' @param leaders Named numeric matrix with columns per leader and rows of their normalized angle disagreement and in-group status.
//' @param dists Named numeric matrix with rows per leader and columns per cell with distances from each cell to chosen cell. 
//' @return Numeric vector of utilities for each cell.
//' @export
// [[Rcpp::export]]
NumericVector flUtility_rcpp(double aFL, double bFL, double dFL, NumericMatrix leaders, NumericMatrix dists) {
  
  // convert without copying
  arma::mat leaders_mat(leaders.begin(), leaders.nrow(), leaders.ncol(), false);
  arma::mat dists_mat(dists.begin(), dists.nrow(), dists.ncol(), false);
  
  // get leader influence
  arma::rowvec b = (bFL + dFL * leaders_mat.row(2)) % leaders_mat.row(1);
  
  // take power of distance matrix
  arma::mat utility = arma::pow(dists_mat, aFL);
  
  // multiply distance by leader influence
  utility.each_col() %= b.t();
  
  // take sum over all leaders
  arma::rowvec col_utility = -sum(utility, 0);
  
  // returns a numeric vector
  return wrap(col_utility);
}


//' Goal Angle Utility
//'
//' @param bGA Numeric scalar scale parameter.
//' @param aGA Numeric scalar power parameter.
//' @param GA Numeric vector of angles to next goal.
//' @param n_rings Integer number of cell rings for which to calculate utlity.
//' 
//' @returns Numeric vector of utilities for each cell.
//' @export
// [[Rcpp::export]]
NumericVector gaUtility_rcpp(double bGA, double aGA, NumericVector GA, int n_rings = 3) {
  
  NumericVector GA_pow = pow(GA, aGA);
  NumericVector output = rep(-(bGA * GA_pow), n_rings);
  
  return output;
}


//' Interpersonal Distance Utility
//' 
//' @param bID Numeric scalar scale parameter.
//' @param dID Numeric scalar direction parameter.
//' @param aID Numeric scalar power parameter.
//' @param is_ingroup Logical vector indicating which subjects in \code{ID_} are part of the ingroup.
//' @param ok Logical matrix indicating if cells are blocked.
//' @param ID_ Numeric matrix or NULL; if not NULL, a numeric matrix of predicted distances from the subject to other pedestrians in the front.
//' @param utility Numeric pre-calculated utility vector (contains \code{-Inf} for blocked cells and \code{0} otherwise).
//' 
//' @returns Numeric vector of utilities for each cell.
//' @export
// [[Rcpp::export]]
NumericVector idUtility_rcpp(double bID, double dID, double aID,
                             LogicalVector is_ingroup,
                             LogicalMatrix ok,
                             Nullable<NumericMatrix> ID_,
                             NumericVector utility) {
  
  
  // None in front return -Inf for cells blocked by objects
  if(is_false(any(ok)) || bID == 0.0) {
    return utility;
  }
    
  NumericMatrix ID(ID_); //re-initialise ID if it is not NULL
  
  NumericVector bID_group = ifelse(is_ingroup, bID, bID + dID);
  
  // convert without copying
  arma::mat ID_mat(ID.begin(), ID.nrow(), ID.ncol(), false);
  // convert with copying and transform into double type
  arma::vec ok_vec = as<arma::vec>(ok);
  // convert without copying
  arma::vec bID_group_vec(bID_group.begin(), bID_group.size(), false);
  // select ok columns and take to power
  arma::mat ID_mat_sub = arma::pow(ID_mat.cols(find(ok_vec)), aID);
  // transform to inverse and muliply column wise
  ID_mat_sub.transform( [](double x) { return (1.0/x); } ).each_col() %= -bID_group_vec;
  // column wise sum
  arma::rowvec col_sums = sum(ID_mat_sub, 0);
  // convert with copying because we don't want to modify utility
  arma::vec utility_vec = as<arma::vec>(utility);
  // replace ok values
  utility_vec.elem(find(ok_vec)) = col_sums;
  
  return wrap(utility_vec);
}


//' Preferred Speed Utility
//'
//' @param aPS Numeric scalar power parameter.
//' @param bPS Numeric scalar scale parameter.
//' @param sPref Numeric scalar preference parameter.
//' @param sSlow Numeric scalar slowness parameter.
//' @param v Numeric scalar current speed.
//' @param d Numeric scalar distance to next goal.
//' @returns Numeric vector of utilities for each cell.
//' @export
// [[Rcpp::export]]
NumericVector psUtility_rcpp(double aPS, double bPS, double sPref, double sSlow,
                             double v, double d) {
  
  // take the parallel min of sPref
  double sPref2 = std::min(sPref, sPref * (d / (v * sSlow)));	
  
  // compute utility for different rings
  NumericVector rings = NumericVector::create(1.5, 1.0, 0.5); // !
  NumericVector output = -bPS * rep_each(pow(abs(v * rings - sPref2), aPS), 11);
  
  return output;
}

//' Walk-beside Utility
//'
//' @param aWB Numeric scalar power parameter.
//' @param bWB Numeric scalar scale parameter.
//' @param buddies Numeric matrix of buddy positions and angles. # needs rewrite
//' @param dists Numeric matrix of distances from cells' centers to closest buddy. # needs rewrite
//' @returns Numeric vector of utilities for each cell.
//' @export
// [[Rcpp::export]]
NumericVector wbUtility_rcpp(double aWB, double bWB, NumericMatrix buddies,
                             NumericMatrix dists) {
  // convert without copying
  arma::mat buddies_mat(buddies.begin(), buddies.nrow(), buddies.ncol(), false);
  arma::mat dists_mat(dists.begin(), dists.nrow(), dists.ncol(), false);
  
  // take power of distance matrix
  arma::mat utility = arma::pow(dists_mat, aWB);
  
  // muliply buddy angle difference with distance
  utility.each_col() %= buddies_mat.row(1).t() * bWB;
  
  // take sum over buddies
  arma::rowvec col_utility = -sum(utility, 0);
  
  // returns a numeric vector
  return wrap(col_utility);
}


//' Total Utility of Cells
//'
//' @param p Numeric vector of subject parameters.
//' @param v Numeric scalar indicating the current speed.
//' @param d Numeric scalar indicating the distance to next goal.
//' @param ba_ NULL or numeric vector of distances from each cell to closest pedestrian.
//' @param ga Numeric vector of angles to next goal.
//' @param id_ NULL or numeric matrix of predicted distances from the subject to other pedestrians in the front.
//' @param is_ingroup Logical vector indicating which subjects in \code{ID_} are part of the ingroup.
//' @param id_pre_utility Numeric pre-calculated interpersonal distance utility vector (contains \code{-Inf} for blocked cells and \code{0} otherwise).
//' @param fl_ NULL or list of numeric matrices:
//' \describe{
//'   \item{leaders}{Numeric matrix of buddy positions and angles.}
//'   \item{dists}{Matrix with rows per leader and columns per cell with distances from each cell to chosen cell.}
//' }
//' @param wb_ NULL or list of numeric matrices:
//' \describe{
//'   \item{buddies}{Matrix with columns per leader and rows of their normalized angle disagreement and in-group status. # needs rewrite}
//'   \item{dists}{Numeric matrix of distances from cells' centers to closest buddy. # needs rewrite}
//' }
//' @param ok Logical matrix indicating if cells are blocked.
//'
//' @return Numeric vector with total utility for each cell.
//' @export
//'
// [[Rcpp::export]]
NumericVector utility(NumericVector p, double v, double d, 
                      Nullable<NumericVector> ba_,
                      NumericVector ga,
                      Nullable<NumericMatrix> id_,
                      LogicalVector is_ingroup,
                      NumericVector id_pre_utility,
                      Nullable<List> fl_,
                      Nullable<List> wb_,
                      LogicalMatrix ok) {
  
  NumericVector ba_utility (33);
  
  if (ba_.isNotNull()) {
    NumericVector ba(ba_);
    IntegerVector ba_names = char2int(ba.names()) - 1; // c++ indexing
    ba_utility = baUtility_rcpp(p["aBA"], p["bBA"], ba, ba_names);
  }
  
  NumericVector ca_utility = caUtility_rcpp(p["aCA"], p["bCA"], p["bCAlr"]);
  
  NumericVector fl_utility (33);
  
  if (fl_.isNotNull()) {
    List fl(fl_);
    fl_utility = flUtility_rcpp(p["aFL"], p["bFL"], p["dFL"], fl["leaders"], fl["dists"]);
  }
  
  NumericVector ga_utility = gaUtility_rcpp(p["bGA"], p["aGA"], ga);

  NumericVector id_utility = id_pre_utility;
  
  if (id_.isNotNull()) {
    id_utility = idUtility_rcpp(p["bID"], p["dID"], p["aID"], is_ingroup, ok, id_, id_pre_utility);
  }

  NumericVector ps_utility = psUtility_rcpp(p["aPS"], p["bPS"], p["sPref"], p["sSlow"], v, d);
  
  NumericVector wb_utility (33);
  
  if (wb_.isNotNull()) {
    List wb(wb_);
    wb_utility = wbUtility_rcpp(p["aWB"], p["bWB"], wb["buddies"], wb["dists"]);
  }
  
  NumericVector total_utility = ps_utility + ga_utility + ca_utility + ba_utility + 
    id_utility + fl_utility + wb_utility;
  
  total_utility.push_front(-p["bS"]);
  
  return total_utility / p["rU"];
}