#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector baUtility_rcpp(double aBA, double bBA, NumericVector BA, IntegerVector idx_BA) {
  NumericVector utility (33);
    
  int k = BA.length();
  
  for(int i = 0; i < k; i++) {
    int idx = idx_BA[i];
    utility[idx] = bBA / pow(std::max(BA[i], 0.0), aBA);
  }
  
  return(utility * -1.0);
}
