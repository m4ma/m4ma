
#ifndef GEOMETRY
#define GEOMETRY

#include <Rcpp.h>

Rcpp::NumericVector dist1_rcpp(Rcpp::NumericVector p1, Rcpp::NumericMatrix p2);

Rcpp::NumericVector angle2_rcpp(Rcpp::NumericMatrix p1, Rcpp::NumericMatrix p2);

Rcpp::NumericMatrix aTOd_rcpp(Rcpp::NumericVector a);

Rcpp::NumericVector Iangle_rcpp(
    Rcpp::NumericMatrix p1,
    double a1,
    Rcpp::NumericMatrix p2
);

Rcpp::NumericVector Dn_rcpp(Rcpp::NumericMatrix p_n, Rcpp::NumericMatrix P_n);

Rcpp::NumericVector minAngle_rcpp(double a1_double, Rcpp::NumericVector a2);

Rcpp::NumericMatrix headingAngle_rcpp(Rcpp::NumericVector a2, double a1);

Rcpp::NumericVector scaleVel_rcpp(Rcpp::NumericVector v, double tStep);

Rcpp::NumericMatrix c_vd_rcpp(
    Rcpp::IntegerVector cells,
    Rcpp::NumericVector p1,
    Rcpp::NumericVector v1,
    double a1, 
    Rcpp::NumericMatrix vels, 
    Rcpp::NumericMatrix angles,
    double tStep
);

Rcpp::NumericMatrix get_vels();

Rcpp::NumericMatrix get_angles();

#endif // GEOMETRY
