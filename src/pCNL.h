
#ifndef PCNL
#define PCNL

#include <Rcpp.h>

double pcnl_rcpp(Rcpp::NumericVector cell, Rcpp::NumericVector utility, 
                 Rcpp::NumericVector mum, Rcpp::List nests, Rcpp::List alpha, 
                 double mu);

double pmnl_rcpp(int cell, Rcpp::NumericVector utility, Rcpp::LogicalMatrix ok);

  
#endif // PCNL
