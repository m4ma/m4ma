
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
    NumericMatrix p1, 
    double a1, 
    double v1, 
    NumericMatrix p2, 
    NumericVector r, 
    List objects
);

Rcpp::Nullable<Rcpp::List> getLeaders_rcpp(
    int n,
    NumericMatrix p_mat,
    NumericVector a,
    NumericVector v,
    NumericMatrix P1,
    NumericVector group,
    NumericMatrix centres,
    List objects,
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