
#ifndef UTILITIES
#define UTILITIES

#include <Rcpp.h>

Rcpp::NumericVector utility(Rcpp::NumericVector p, int n, double v, double d, 
                            Rcpp::Nullable<Rcpp::NumericVector> ba_,
                            Rcpp::NumericVector ga,
                            Rcpp::Nullable<Rcpp::NumericMatrix> id_,
                            Rcpp::Nullable<Rcpp::List> fl_,
                            Rcpp::Nullable<Rcpp::List> wb_,
                            Rcpp::LogicalMatrix ok,
                            Rcpp::IntegerVector group);

#endif // UTILITIES
