#ifndef UTILS
#define UTILS

#include <Rcpp.h>

Rcpp::IntegerVector char2int(Rcpp::CharacterVector char_vec);

Rcpp::CharacterVector int2char(Rcpp::IntegerVector int_vec);

Rcpp::NumericVector bin_vector(Rcpp::NumericVector x, Rcpp::NumericVector bins);

#endif // UTILS