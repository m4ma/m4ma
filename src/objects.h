
#ifndef OBJECTS
#define OBJECTS

#include <Rcpp.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::d2::point_xy<double> point_t;
typedef bg::model::box<point_t> box_t;
typedef bg::model::polygon<point_t> polygon_t;
typedef bg::model::multi_polygon<polygon_t> multi_polygon_t;
typedef std::pair<box_t, unsigned> rtree_elem_t;
typedef bgi::rtree<rtree_elem_t, bgi::quadratic<16>> rtree_t;

template <> point_t Rcpp::as(SEXP ptsexp);
template <> multi_polygon_t Rcpp::as(SEXP osexp);

rtree_t objects_to_rtree(Rcpp::List objects);

#endif // OBJECTS