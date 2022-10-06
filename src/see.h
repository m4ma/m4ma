
#ifndef SEE
#define SEE

#include <Rcpp.h>

Rcpp::LogicalVector seesMany_rcpp(
    Rcpp::NumericVector p1,
    Rcpp::NumericMatrix ps,
    Rcpp::List objects
);

Rcpp::NumericVector line_line_intersection_rcpp(
    Rcpp::NumericVector P1,
    Rcpp::NumericVector P2,
    Rcpp::NumericVector P3,
    Rcpp::NumericVector P4,
    bool interior_only=false
);

bool seesGoal_rcpp(
    Rcpp::NumericVector p_n,
    Rcpp::NumericVector P_n,
    Rcpp::List objects
);

#endif // SEE