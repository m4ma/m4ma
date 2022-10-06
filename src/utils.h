#ifndef UTILS
#define UTILS

#include <Rcpp.h>

Rcpp::IntegerVector char2int(Rcpp::CharacterVector char_vec);

Rcpp::CharacterVector int2char(Rcpp::IntegerVector int_vec);
  
#endif // UTILS