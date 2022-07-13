#include <Rcpp.h>
using namespace Rcpp;
//' dist_rcpp
//'
//' Compute the istance from p1 to p2, both xy column matrices
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


