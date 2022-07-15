#include <Rcpp.h>
#include <math.h>

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
  int n_rows = p1.nrow();
  NumericVector angle(n_rows);

  for(int row = 0; row < n_rows; ++row){
    double angle_centered = (180.0 / M_PI) * std::atan2((p2(_, 1) - p1(_, 1))[row], (p2(_, 0) - p1(_, 0))[row]);
    angle[row] = fmod(360.0 + angle_centered, 360.0);
    // if(angle_centered < 0.0) {
    //   angle[row] = 360.0 + angle_centered;
    // } else {
    //   angle[row] = angle_centered;
    // }
  }

  return angle;
}
