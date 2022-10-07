#include <Rcpp.h>

using namespace Rcpp;

IntegerVector char2int(CharacterVector char_vec) {
  int k = char_vec.length();
  
  IntegerVector int_vec (k);
  
  for(int i = 0;  i < k; ++i) {
    int_vec[i] = std::stoi(std::string(char_vec[i]));
  }
  
  return int_vec;
}


CharacterVector int2char(IntegerVector int_vec) {
  int k = int_vec.length();
  
  CharacterVector char_vec (k);
  
  for(int i = 0;  i < k; ++i) {
    char_vec[i] = std::to_string(int_vec[i]);
  }
  
  return char_vec;
}


// Helper function that bins vector according to borders and returns bin indices
NumericVector bin_vector(NumericVector x, NumericVector bins) {
  NumericVector out(x.length());
  
  for(int i = 0; i < x.length(); i++) {
    // set to NA if outside borders
    double idx = NA_REAL;
    
    // iterate over borders
    for(int j = 1; j < bins.length(); j++) {
      // if angle is > left border and <= right border
      if(x[i] > bins[j-1] && x[i] <= bins[j]) {
        idx = j;
      }
      out[i] = idx;
    }
  }
  
  return(out);
}