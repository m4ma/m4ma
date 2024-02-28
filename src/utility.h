
#ifndef UTILITY
#define UTILITY

#include <Rcpp.h>

Rcpp::NumericVector utility(Rcpp::NumericVector p, double v, double d, 
                            Rcpp::Nullable<Rcpp::NumericVector> ba_,
                            Rcpp::NumericVector ga,
                            Rcpp::Nullable<Rcpp::NumericMatrix> id_,
                            Rcpp::LogicalVector is_ingroup,
                            Rcpp::NumericVector id_pre_utility,
                            Rcpp::Nullable<Rcpp::List> fl_,
                            Rcpp::Nullable<Rcpp::List> wb_,
                            Rcpp::LogicalMatrix ok);

#endif // UTILITY
