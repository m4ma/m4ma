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