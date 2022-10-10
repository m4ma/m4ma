
#ifndef GEOMETRY
#define GEOMETRY

#include <Rcpp.h>

Rcpp::NumericVector dist1_rcpp(Rcpp::NumericVector p1, Rcpp::NumericMatrix p2);

Rcpp::NumericVector angle2_rcpp(Rcpp::NumericMatrix p1, Rcpp::NumericMatrix p2);

Rcpp::NumericMatrix aTOd_rcpp(Rcpp::NumericVector a);

Rcpp::NumericVector minAngle_rcpp(double a1_double, Rcpp::NumericVector a2);

Rcpp::NumericMatrix headingAngle_rcpp(Rcpp::NumericVector a2, double a1);

#endif // GEOMETRY
