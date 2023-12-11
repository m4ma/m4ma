
#ifndef OBJECTS
#define OBJECTS

#include <Rcpp.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/geometries.hpp>

namespace bg = boost::geometry;

typedef bg::model::d2::point_xy<double> point_t;
typedef bg::model::polygon<point_t> polygon_t;
typedef bg::model::multi_polygon<polygon_t> multi_polygon_t;

template <> point_t Rcpp::as(SEXP ptsexp);
template <> multi_polygon_t Rcpp::as(SEXP osexp);

#endif // OBJECTS