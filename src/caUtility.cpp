#include <Rcpp.h>
using namespace Rcpp;

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
