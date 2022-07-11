#include <Rcpp.h>
using namespace Rcpp;
//' Compute flUtility
//'
//' @param aFL numeric vector
//' @param bFL numeric vector
//' @param dFL numeric vector
//' @param leaders named numeric matrix 
//' @param dists transposed numeric matrix 
//' @return numeric vector of 33 elements
//' @examples
//' flUtility_rcpp(1.5, 1, 4, leaders = matrix(c(9, 0.5, 0.5, 8, 0.5, 0.5), ncol = 2), dists = t(matrix(rep(runif(33), 2), ncol = 2)))
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
    col_utility[j] = sum(utility(_, j));
  }
  
  return col_utility * -1.0;
}
