#include <Rcpp.h>
using namespace Rcpp;
//' dist1_rcpp
//'
//' Compute distance from p1 to p2
//' @param p1 Numeric matrix of a single point (i.e., 1 row, 2 xy columns)
//' @param p2 Numeric matrix of multiple xy points
//' @return Named numeric vector of length equal to the number of rows in p2
//' @export
// [[Rcpp::export]]
NumericVector dist1_rcpp(NumericMatrix p1, NumericMatrix p2) {
  
  // allocate transposed matrix
  NumericMatrix t_p2 = transpose(p2);
  // allocate indices
  int ncols = t_p2.ncol();
  int nrows = t_p2.nrow();
  // allocate output
  NumericMatrix temp_matrix(t_p2);
  NumericVector out(ncols); 
  
  for(int i = 0; i < nrows; i++){
    for(int j = 0; j < ncols; j++){
      temp_matrix(i,j) = std::pow(t_p2(i,j) - p1(0,i),2);
      out(j) = std::sqrt(sum(temp_matrix(_, j)));
    }
  }
  return out;
}
