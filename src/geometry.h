
#ifndef GEOMETRY
#define GEOMETRY

#include <Rcpp.h>

Rcpp::NumericVector geometry(Rcpp::NumericMatrix p1, Rcpp::NumericMatrix p2, 
                             Rcpp::NumericVector a, Rcpp::NumericVector a1,
                             double a1_double, double tStep, double v,
                             Rcpp::NumericVector a2, Rcpp::NumericMatrix vels, 
                             Rcpp::NumericMatrix angles);

#endif // GEOMETRY
