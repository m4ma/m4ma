#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

//' angle2
//' 
//' Anti-clockwise angle from p1 as origin to p2 (x,y pairs matrices) 0 to <360
//' @param p1 Numeric matrix
//' @param p2 Numeric matrix
//' @return Named numeric vector of length equal to the number of rows in p1


// [[Rcpp::export]]
NumericVector angle2_rcpp(NumericMatrix p1, NumericMatrix p2) {
  // initialise variable
  int nrow = p1.nrow();
  NumericVector out(nrow);
  NumericVector temp(nrow);
  
  for(int row = 0; row < nrow; row++){
    temp(row) = (180 / M_PI) * std::atan2((p2(_,1) - p1(_,1))[row], (p2(_,0) - p1(_,0))[row]) ;
    Rcout << "The value of temp : " << temp(row) << "\n";
    out(row) = std::fmodl(temp(row), 360);
    Rcout << "The value of out : " << out << "\n";
    
  }
  
  return(out);
}
