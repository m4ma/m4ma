// [[Rcpp::depends(BH)]]

// Enable C++14 via this plugin to compile boost
// [[Rcpp::plugins("cpp14")]]

#include <Rcpp.h>
#include <math.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/index/rtree.hpp>
// #include <boost/geometry/geometries/adapted/boost_polygon.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
// namespace trans = boost::geometry::strategy::transform;

namespace Rcpp {

typedef bg::model::d2::point_xy<double> point_t;
typedef bg::model::box<point_t> box_t;
typedef bg::model::polygon<point_t> polygon_t;
typedef bg::model::multi_polygon<polygon_t> multi_polygon_t;
typedef std::pair<box_t, unsigned> rtree_elem_t;
typedef bgi::rtree<rtree_elem_t, bgi::rstar<16>> rtree_t;
// typedef boost::polygon::rectangle_data<double> rect;

template <> point_t as(SEXP ptsexp) {
  Rcpp::NumericVector pt(ptsexp);
  return point_t(pt[0], pt[1]);
};

template <> SEXP wrap(const point_t &p) {
  return wrap(Rcpp::NumericVector::create(p.x(), p.y()));
};

// Template for rectangle objects
// template <> rect as(SEXP rsexp) {
//   Rcpp::S4 b(rsexp);
//   
//   NumericVector p1 = b.slot("lower");
//   NumericVector p2 = b.slot("upper");
//   
//   rect r = boost::polygon::construct<rect>(p1[0], p1[1], p2[0], p2[1]);
//   
//   rect r_t;
//   
//   double angle = b.slot("orientation");
//   
//   trans::rotate_transformer<boost::geometry::radian, double, 2, 2> rotate(-angle);
//   boost::geometry::transform(r, r_t, rotate);
//   
//   return r_t;
// };

multi_polygon_t s4_circle_to_multi_polygon(S4 obj) {
  Rcpp::NumericVector center = obj.slot("center");
  double radius = obj.slot("radius");
  
  point_t p = as<point_t>(center);
  
  multi_polygon_t circle;
  
  const int n_points = 36;
  
  bg::strategy::buffer::point_circle point_strategy(n_points);
  bg::strategy::buffer::distance_symmetric<double> distance_strategy(radius);
  bg::strategy::buffer::join_round join_strategy(n_points);
  bg::strategy::buffer::end_round end_strategy(n_points);
  bg::strategy::buffer::side_straight side_strategy;
  
  bg::buffer(p, circle, distance_strategy, side_strategy, join_strategy, end_strategy, point_strategy);
  
  return circle;
}

multi_polygon_t s4_polygon_to_multi_polygon(S4 obj) {
  Rcpp::NumericMatrix points = obj.slot("points");
  
  multi_polygon_t mpoly;
  
  mpoly.resize(1);
  
  for (int i = 0; i < points.nrow(); i++) {
    Rcpp::NumericVector c_i = points(i, _);
    point_t p_i = point_t(c_i[0], c_i[1]);
    
    bg::append(mpoly[0].outer(), p_i);
  }
  
  return mpoly;
}

multi_polygon_t list_to_multi_polygon(List obj) {
  Rcpp::NumericVector x = obj["x"];
  Rcpp::NumericVector y = obj["y"];
  // create object box from min_corner, max_corner
  box_t bx(
      point_t(x[0], y[0]),
      point_t(x[1], y[1])
  );
  
  multi_polygon_t mpoly;
  
  bg::convert(bx, mpoly);
  
  return mpoly;
}

template <> multi_polygon_t as(SEXP osexp) {
  Rcpp::RObject r_obj(osexp);
  
  multi_polygon_t mpoly;
  
  if (r_obj.isS4() && r_obj.inherits("circle")) {
    Rcpp::S4 obj(r_obj);
    mpoly = s4_circle_to_multi_polygon(obj);
  } else if (r_obj.isS4() && r_obj.inherits("polygon")) {
    Rcpp::S4 obj(r_obj);
    mpoly = s4_polygon_to_multi_polygon(obj);
  } else if (is<List>(r_obj)) {
    Rcpp::List obj(r_obj);
    mpoly = list_to_multi_polygon(obj);
  }
  
  return mpoly;
};

}

using namespace Rcpp;

rtree_t objects_to_rtree(List objects) {
  std::vector<rtree_elem_t> boxes;
  
  for (unsigned i = 0; i < objects.length(); i++) {
    RObject objects_i = objects[i];
    
    // instantiate object as a polygon
    multi_polygon_t mpoly_i = as<multi_polygon_t>(objects_i);
    
    // rtree requires bounding box as input
    box_t box_i = bg::return_envelope<box_t>(mpoly_i);
    
    boxes.push_back(std::make_pair(box_i, i));
  }
  
  // construct rtree
  rtree_t rtree(boxes);
  
  return rtree;
}
