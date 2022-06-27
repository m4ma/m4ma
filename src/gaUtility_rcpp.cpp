#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]


NumericVector gaUtility_rcpp(int bGA, int aGA, NumericVector GA) {
  
  NumericVector GA_pow = pow(GA, aGA);
  NumericVector output = rep((bGA * GA_pow) * -1, 3);
  
  return output;
}