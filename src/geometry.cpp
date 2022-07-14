#include <Rcpp.h>

using namespace Rcpp;

//' dist_rcpp
//'
//' Compute distance from p1 to p2, both xy column matrices
//' @param p1 Numeric matrix
//' @param p2 Numeric matrix
//' @return Named numeric vector of length equal to the number of rows in p1
//' @export
// [[Rcpp::export]]
NumericVector dist_rcpp(NumericMatrix p1, NumericMatrix p2) {
  
  // initialise variables
  int ncols = p1.ncol();
  int nrows = p1.nrow();
  NumericMatrix temp(p1);
  CharacterVector p1_names = rownames(p1);
  NumericVector out(nrows);
  out.names() = p1_names;
  
  // difference between the rows of p1 and p2 elevated to the power of 2
  for(int i = 0; i < nrows; i++){
    for(int j = 0; j < ncols; j++){
      temp(i,j) = std::pow(p1(i,j) - p2(i,j),2);
    }
  }
  
  // apply square root of the sum of the rows of x and y
  for(int i = 0; i < nrows; i++){
    out(i) = std::sqrt(temp(i,0) + temp(i,1));
  }
  
  return out;
}


//' dist1_rcpp
//'
//' Compute distance from p1 to p2
//' @param p1 Numeric matrix of a single point (i.e., 1 row, 2 xy columns)
//' @param p2 Numeric matrix of multiple xy points
//' @return Named numeric vector of length equal to the number of rows in p2
//' @export
// [[Rcpp::export]]
NumericVector dist1_rcpp(NumericMatrix p1, NumericMatrix p2) {
  
  // allocate transposed matrix
  NumericMatrix t_p2 = transpose(p2);
  // allocate indices
  int ncols = t_p2.ncol();
  int nrows = t_p2.nrow();
  // allocate temporary matrix
  NumericMatrix temp_matrix(t_p2);
  // allocate output
  NumericVector out(ncols); 
  
  for(int i = 0; i < nrows; i++){
    for(int j = 0; j < ncols; j++){
      temp_matrix(i,j) = std::pow(t_p2(i,j) - p1(0,i),2);
      out(j) = std::sqrt(sum(temp_matrix(_, j)));
    }
  }
  return out;
}


//' angle2s_rcpp
//'
//' Compute Shortest angle anti-clockwise from p1 as origin to p2 (> -180 to 180)
//' @param p1 Numeric matrix
//' @param p2 Numeric matrix
//' @return Numeric vector of length equal to the number of rows in p1 
//' @export
// [[Rcpp::export]]
NumericVector angle2s_rcpp(NumericMatrix p1, NumericMatrix p2) {
  
  // initialise variable
  NumericVector out(p1.nrow());
  CharacterVector p1_names = rownames(p1);
  int nrow = p1.nrow();
  
  for(int row = 0; row < nrow; row++){
    out(row) = (180 / M_PI) * std::atan2((p2(_,1) - p1(_,1))[row], (p2(_,0) - p1(_,0))[row]);
  }
  
  out = round(out, 10);
  
  // return a numeric vector whose names match those of p1's rows
  // only if p1's rows have names
  if(is_false(all(is_na(p1_names)))){
    out.names() = p1_names;
  }
  
  return(out);
}
