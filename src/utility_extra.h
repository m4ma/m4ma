
#ifndef UTILITY_EXTRA
#define UTILITY_EXTRA

#include <Rcpp.h>

Rcpp::NumericVector destinationAngle_rcpp(
    double a,
    Rcpp::NumericMatrix p1,
    Rcpp::NumericMatrix P1
);

Rcpp::Nullable<Rcpp::NumericMatrix> predClose_rcpp(
    int n,
    Rcpp::NumericMatrix p1,
    double a1,
    Rcpp::NumericMatrix p2,
    Rcpp::NumericVector r,
    Rcpp::NumericMatrix centres,
    Rcpp::NumericMatrix p_pred,
    Rcpp::List objects
);

Rcpp::NumericVector blockedAngle_rcpp(
    int n,
    Rcpp::List state,
    Rcpp:: NumericMatrix p_pred,
    Rcpp::List objects
);

Rcpp::Nullable<Rcpp::List> getLeaders_rcpp(
    int n,
    Rcpp::List state,
    Rcpp::NumericMatrix centres,
    Rcpp::List objects,
    bool onlyGroup = false,
    bool preferGroup = true,
    bool pickBest = false
);

Rcpp::Nullable<Rcpp::List> getBuddy_rcpp(
    int n,
    Rcpp::NumericVector group,
    Rcpp::NumericVector a,
    Rcpp::NumericMatrix p_pred,
    Rcpp::NumericMatrix centres,
    Rcpp::List objects,
    bool pickBest,
    Rcpp::List state
);

#endif // UTILITY_EXTRA