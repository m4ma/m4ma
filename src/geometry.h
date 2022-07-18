
#ifndef GEOMETRY
#define GEOMETRY

#include <Rcpp.h>

Rcpp::NumericVector geometry(Rcpp::NumericMatrix p1, Rcpp::NumericMatrix p2, 
                             Rcpp::NumericVector a, Rcpp::NumericVector a1,
                             double a1_double, double tStep,
                             Rcpp::NumericVector a2, Rcpp::NumericVector v);

#endif // GEOMETRY
