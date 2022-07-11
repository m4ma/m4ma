#include <Rcpp.h>
using namespace Rcpp;

//' Compute gaUtility (Goal angle utility)
//'
//' @param bGA integer
//' @param aGA integer
//' @param GA numeric vector
//' @returns a numeric vector of length equal to GA's length 
//' @examples
//' gaUtility_rcpp(bGA = 10, aGA = 2, GA = rnorm(11, mean = .54, sd = 0.38))
//' @export
// [[Rcpp::export]]

NumericVector gaUtility_rcpp(int bGA, int aGA, NumericVector GA) {

  NumericVector output(3 * GA.length());
  
  for (int ix = 0; ix < GA.length(); ++ix) {
    auto number = -bGA * std::pow(GA[ix], aGA);
    output[ix] = number;
    output[ix + GA.length()] = number;
    output[ix + 2 * GA.length()] = number;
  }

  return output;
}