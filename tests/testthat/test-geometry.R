# Test rcpp geometry functions against R counterparts

# Geometry functions take as inputs angles, or xy coordinates

set.seed(123)

p2 = matrix(rnorm(6), ncol = 2, nrow = 3)
rownames(p2) = c('a_i', 'b_i', 'c_i')
colnames(p2) = c('x', 'y')

p1 = matrix(rnorm(6), ncol = 2, nrow = 3)
rownames(p1) = c('a_i', 'b_i', 'c_i')
colnames(p1) = c('x', 'y')

single_point = matrix(rnorm(2), ncol = 2, nrow = 1)
rownames(single_point) = c('a_i')
colnames(single_point) = c('x', 'y')

a = c(100, 23, 300)
names(a) = c('a_i', 'b_i', 'c_i')

a_double = c(a_i = -90)

v = c(a_i = 1)

tStep = 0.5


testthat::test_that("geometry wrapper works", {
  fun = m4ma::wrapper('dist1')
  res_rcpp = fun(single_point, p2, use = 'cpp')
  res_r = fun(single_point, p2, use = 'r')
  testthat::expect_equal(res_rcpp, res_r)
})


testthat::test_that("dist geometry function works", {
  ref = m4ma::dist_r(p1, p2)
  res = m4ma::dist_rcpp(p1, p2)
  testthat::expect_equal(res, ref)
})

testthat::test_that("dist1 geometry function works", {
  ref = m4ma::dist1_r(single_point, p2)
  res = m4ma::dist1_rcpp(single_point, p2)
  testthat::expect_equal(res, ref)
})

testthat::test_that("angle2 geometry function works (single point)", {
  ref = m4ma::angle2_r(single_point, p2)
  res = m4ma::angle2_rcpp(single_point, p2)
  testthat::expect_equal(res, ref)
})

testthat::test_that("angle2 geometry function works (matrix)", {
  ref = m4ma::angle2_r(p1, p2)
  res = m4ma::angle2_rcpp(p1, p2)
  testthat::expect_equal(res, ref)
})

testthat::test_that("angle2s geometry function works (single point)", {
  ref = m4ma::angle2s_r(single_point, p2)
  res = m4ma::angle2s_rcpp(single_point, p2)
  testthat::expect_equal(res, ref)
})

testthat::test_that("angle2s geometry function works (matrix)", {
  ref = m4ma::angle2s_r(p1, p2)
  res = m4ma::angle2s_rcpp(p1, p2)
  testthat::expect_equal(res, ref)
})

testthat::test_that("aTOd geometry function works", {
  ref = m4ma::aTOd_r(a)
  res = m4ma::aTOd_rcpp(a)
  testthat::expect_equal(res, ref)
})

testthat::test_that("Iangle geometry function works", {
  ref = m4ma::Iangle_r(single_point, a_double, p2)
  res = m4ma::Iangle_rcpp(single_point, a_double, p2)
  testthat::expect_equal(res, ref)
})

testthat::test_that("Dn geometry function works", {
  ref = m4ma::Dn_r(p_n = p1, P_n = p2) 
  res = m4ma::Dn_rcpp(p_n = p1, P_n = p2)
  testthat::expect_equal(unname(res), ref)
})

testthat::test_that("minAngle geometry function works", {
  ref = m4ma::minAngle_r(a_double, a) 
  res = m4ma::minAngle_rcpp(a_double, a)
  testthat::expect_equal(res, ref)
})

testthat::test_that("scaleVel geometry function works", {
  ref = m4ma::scaleVel_r(v, tStep) 
  res = m4ma::scaleVel_rcpp(v, tStep)
  testthat::expect_equal(res, ref)
})

testthat::test_that("coneNum geometry function works", {
  ref = m4ma::coneNum_r(1:33) 
  res = m4ma::coneNum_rcpp(1:33)
  testthat::expect_equal(res, ref)
})

testthat::test_that("ringNum geometry function works", {
  ref = m4ma::ringNum_r(1:33) 
  res = m4ma::ringNum_rcpp(1:33)
  testthat::expect_equal(res, ref)
})

vels = matrix(rep(c(1.5, 1, .5), each = 11), ncol = 3)
angles = matrix(rep(c(72.5, 50, 32.5, 20, 10, 0, 350, 340, 
                      327.5, 310, 287.5), times = 3), ncol = 3)

testthat::test_that("c_vd geometry function works", {
  ref = m4ma::c_vd_r(cells = 1:33, p1 = single_point[1, ],
                   v1 = v, a1 = a_double) 
  res = m4ma::c_vd_rcpp(cells = 1:33, p1 = single_point[1, ],
                        v1 = v, a1 = a_double,
                        vels = vels,
                        angles = angles,
                        tStep = 0.5)
  testthat::expect_equal(res, ref)
})




