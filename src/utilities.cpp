#include <Rcpp.h>

using namespace Rcpp;


//' Blocked-angle Utility
//'
//' @param aBA Power parameter.
//' @param bBA Scale parameter.
//' @param BA Numeric vector of distances from cells to closest pedestrians.
//' @param idx_BA Integer vector of cell indices.
//'
//' @return Numeric Vector of utilities for each cell.
//' @export
// [[Rcpp::export]]
NumericVector baUtility_rcpp(double aBA, double bBA, NumericVector BA, IntegerVector idx_BA) {
  NumericVector utility (33);
  
  int k = BA.length();
  
  for(int i = 0; i < k; i++) {
    int idx = idx_BA[i];
    utility[idx] = -bBA / std::pow(std::max(BA[i], 0.0), aBA);
  }
  
  return(utility);
}


//' Compute caUtility
//' @param aCA numeric vector
//' @param bCA numeric vector
//' @param bCAlr numeric vector
// [[Rcpp::export]]
NumericVector caUtility_rcpp(double aCA, double bCA, double bCAlr) {
  
  // allocate variables
  NumericVector output(11);
  
  // angles 
  NumericVector angles =
    NumericVector::create(10, 20, 32.5, 50, 72.5)/90 ;
  
  // compute power of angles
  NumericVector ap = pow(angles, aCA);
  
  // fill the output
  int ap_len = ap.length();
  
  for(int i = 0;  i < ap_len; ++i) {
    output[i] = bCA * bCAlr *ap[ap_len-1-i];
    output[i+ap_len+1] = (bCA / bCAlr) *ap[i];
  }
  output[ap_len] = 0;
  
  
  output = -rep(output, 3);
  return output;
  
}


//' Compute flUtility
//'
//' @param aFL numeric vector
//' @param bFL numeric vector
//' @param dFL numeric vector
//' @param leaders named numeric matrix 
//' @param dists transposed numeric matrix 
//' @return numeric vector of 33 elements
//' @export
// [[Rcpp::export]]

NumericVector flUtility_rcpp(double aFL, double bFL, double dFL, NumericMatrix leaders, NumericMatrix dists) {
  
  NumericVector b = (bFL + dFL * leaders(2, _)) * leaders(1, _);
  
  int k = b.length();
  
  NumericMatrix utility (k, dists.ncol());
  
  for(int i = 0; i < k; ++i) {
    utility(i, _) = b[i] * pow(dists(i, _), aFL);
  }
  
  int cols = utility.ncol();
  
  NumericVector col_utility (cols);
  
  for(int j = 0; j < cols; ++j) {
    col_utility[j] = -sum(utility(_, j));
  }
  
  return col_utility;
}


//' Compute gaUtility (Goal angle utility)
//'
//' @param bGA integer
//' @param aGA integer
//' @param GA numeric vector
//' @returns a numeric vector of length equal to GA's length 
//' @export
// [[Rcpp::export]]
NumericVector gaUtility_rcpp(double bGA, double aGA, NumericVector GA) {
  
  NumericVector GA_pow = pow(GA, aGA);
  NumericVector output = rep(-(bGA * GA_pow), 3);
  
  return output;
}


//' Compute IDUtility
//'
//' Inter-personal distance utility for cell 1..33. b parameter divided by
//' sum over power of distances between bodies for cell to all inFront peds.
//' 
//' @param bID Numeric scalar.
//' @param dID Numeric scalar.
//' @param aID Numeric scalar.
//' @param n Numeric scalar indexing the subject in the state.
//' @param ok Logical matrix indicating if cells are blocked.
//' @param group Named numeric scalar with group indices for each pedestrian.
//' @param ID_ Numeric matrix of the type NULL - if not NULL, ID_ is a Numeric matrix of predicted distances from the subject to other pedestrians in the front.
//' @returns Numeric vector with interpersonal distance utilities for each cell.
//' @export
// [[Rcpp::export]]
NumericVector idUtility_rcpp(double bID, double dID, double aID, double n, 
                             const LogicalMatrix ok, IntegerVector group,
                             Nullable<NumericMatrix> ID_) {
  
  // allocate output
  int ok_rows = ok.nrow();
  int ok_cols = ok.ncol();
  NumericVector utility(ok_rows * ok_cols);
  
  // None in front return -Inf for cells blocked by objects
  if(Rf_isNull(ID_) || is_false(any(ok))) {
    for(int row = 0; row < ok_rows; row++) {
      for(int col = 0; col < ok_cols; col++) {
        if(ok(row, col)){
          utility[col*ok_rows+row] = 0;
        } else{
          utility[col*ok_rows+row] = R_NegInf;
        }
      }
    }
    
    return utility;
  }
    
  NumericMatrix ID(ID_); //re-initialise ID if it is not NULL
  int ID_rows = ID.nrow();
  int ID_cols = ID.ncol();
  
  CharacterVector ID_row_names = rownames(ID);
  
  IntegerVector group_excl_n = group;
  group_excl_n.erase(n); // remove nt element
  
  IntegerVector in_group = group_excl_n[group_excl_n == group[n]];
  
  CharacterVector names_in_group = {""};
  
  if (in_group.isNULL()) {
    
  } else {
    CharacterVector names_in_group = in_group.names();
  }
  NumericVector bID_2 = (ID_rows);
  
  // Group dependent b, bigger for outgroup by dID
  for(int row = 0; row < ID_rows; row++) {
    if(names_in_group.containsElementNamed(ID_row_names[row])) {
      
      bID_2[row] = bID;
      
    } else {
      
      bID_2[row] = bID + dID;
    }
  }
  
  // Object or pedestrian clash
  LogicalVector ID_cols_gt_0 (ID_cols);
  
  for(int col = 0; col < ID_cols; col++) {
    ID_cols_gt_0[col] = is_true(all(ID(_, col) > 0.0));
  }
  
  LogicalMatrix ok_clash (ok);
  
  for(int row = 0; row < ok_rows; row++) {
    for(int col = 0; col < ok_cols; col++) {
      ok_clash(row, col) = ok(row, col) & ID_cols_gt_0[col*ok_rows+row]; // column-wise fill
      if(ok_clash(row, col)) {
        utility[col*ok_rows+row] = 0;
      } else {
        utility[col*ok_rows+row] = R_NegInf;
      }
    }
  }
  
  // Repulsion
  if(bID != 0) {
    for(int col = 0; col < ID_cols; col++) {
      if(ok_clash[col]) {
        NumericVector col_utility (ID_rows); 
        for(int row = 0; row < ID_rows; row++) {
          col_utility[row] = -(bID_2[row] / (std::pow(ID(row, col), aID)));
        }
        utility[col] = sum(col_utility);
      }
    }
  }
  
  return utility;
}


//' Compute psUtility
//'
//' @param aPS numeric vector
//' @param sPref numeric vector
//' @param sSlow numeric vector
//' @param bPS numeric vector
//' @param v numeric vector
//' @param d numeric vector 
//' @returns a numeric vector of length equals to d's length 
//' @export
// [[Rcpp::export]]
NumericVector psUtility_rcpp(double aPS, double bPS, double sPref, double sSlow,
                             double v, double d) {
  // allocate variables
  NumericVector output(33);
  
  // take the parallel min of sPref
  double sPref2 = std::min(sPref, d * sPref / sSlow);	
  
  // fill output vector 
  for(int i = 0;  i < 11; ++i) {
      output[i] = -bPS * abs(pow(v * 1.5 - sPref2, aPS));
      output[(i+11)] = -bPS * abs(pow(v - sPref2, aPS));
      output[(i+22)] = -bPS * abs(pow(v/2 - sPref2, aPS));  
  }
  
  return output;
  
}

// [[Rcpp::export]]
NumericVector wbUtility_rcpp(double aWB, double bWB, NumericMatrix buddies,
                             NumericMatrix dists) {
  
  int n_rows = dists.rows();
  int n_cols = dists.cols();
  NumericVector utility (n_cols);
  
  for(int i = 0;  i < n_cols; ++i) {
    NumericVector col_utility (n_rows);
    for(int j = 0;  i < n_rows; ++i) {
      col_utility = bWB * buddies(1, j) * pow(dists(j, _), aWB);
    }
    utility[i] = -sum(col_utility);
  }
  
  return utility;
}


IntegerVector char2int(CharacterVector char_vec) {
  int k = char_vec.length();
  
  IntegerVector int_vec (k);
  
  for(int i = 0;  i < k; ++i) {
    int_vec[i] = std::stoi(std::string(char_vec[i]));
  }
  
  return int_vec;
}


//' Total Utility of Cells
//'
//' @param p Numeric vector of subject parameters.
//' @param n Integer scalar indexing the subject in the state.
//' @param v Numeric scalar indicating the current speed.
//' @param d Numeric scalar indicating the distance to next goal.
//' @param ba_ Numeric vector of distances from each cell to closest pedestrian.
//' @param ga Numeric vector of angles to next goal.
//' @param id_ Numeric matrix of predicted distances from the subject to other pedestrians in the front.
//' @param fl_ List of numeric matrices:
//' \describe{
//'   \item{leaders}{Matrix with columns per leader and rows of their normalized angle disagreement and in-group status.}
//'   \item{dists}{Matrix with rows per leader and columns per cell with distances from each cell to chosen cell.}
//' }
//' @param wb_ Numeric vector of distances from cells' centers to closest buddy.
//' @param ok Logical matrix indicating if cells are blocked.
//' @param group Integer vector with group indices for each pedestrian.
//'
//' @return Numeric vector with total utility for each cell.
//' @export
//'
// [[Rcpp::export]]
NumericVector utility(NumericVector p, int n, double v, double d, 
                           Nullable<NumericVector> ba_,
                           NumericVector ga,
                           Nullable<NumericMatrix> id_,
                           Nullable<List> fl_,
                           Nullable<List> wb_,
                           LogicalMatrix ok,
                           IntegerVector group) {
  
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

  NumericVector id_utility = idUtility_rcpp(p["bID"], p["dID"], p["aID"], n, ok, group, id_);

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