#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

//' angle2s_rcpp
//'
//' Compute Shortest angle anti-clockwise from p1 as origin to p2 (> -180 to 180)
//' @param p1 Numeric matrix
//' @param p2 Numeric matrix
//' @return Named numeric vector of length equal to the number of rows in p1
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
  if(is_true(all(is_na(p1_names)))){
    return(out);
  } else {
    out.names() = p1_names;
    return(out);
  }
  
  
}


