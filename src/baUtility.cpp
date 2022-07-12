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
    int idx = idx_BA[i]-1;
    utility[idx] = -bBA / std::pow(std::max(BA[i], 0.0), aBA);
  }
  
  return(utility);
}
