#include <Rcpp.h>
using namespace Rcpp;

//' Compute gaUtility (Goal angle utility)
//'
//' @param bGA integer
//' @param aGA integer
//' @param GA numeric vector
//' @returns a numeric vector of length equal to GA's length 
//' @export
// [[Rcpp::export]]

NumericVector gaUtility_rcpp(int bGA, int aGA, NumericVector GA) {

  NumericVector GA_pow = pow(GA, aGA);
  NumericVector output = rep(-(bGA * GA_pow), 3);

  return output;
}