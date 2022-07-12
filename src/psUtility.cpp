#include <Rcpp.h>

using namespace Rcpp;

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

NumericVector psUtility_rcpp(double aPS, double sPref, 
                             double sSlow, double bPS, 
                             double v, NumericVector d) {
  // allocate variables
  int I = d.length();
  NumericVector sPref2(I);
  NumericVector output(33*I);
  
  // take the parallel min of sPref
  for(int x=0; x<I; ++x){
    sPref2[x] = std::min(sPref, d[x] * sPref / sSlow);	
  }
  
  // fill output vector 
  for(int i = 0;  i < 11; ++i) {
    for (int j = 0;  j < I; ++j) {
      
      output[i*I+j] = -bPS * abs(pow(v * 1.5 - sPref2[j], aPS));
      output[(i+11)*I+j] = -bPS * abs(pow(v - sPref2[j], aPS));
      output[(i+22)*I+j] = -bPS * abs(pow(v/2 - sPref2[j], aPS));  
    }
  }
  
  return output;
  
}

