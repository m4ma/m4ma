# Test rcpp geometry functions against R counterparts

# geometry functions take as inputs angles, or xy coordinates

set.seed(123)
p2 <- matrix(rnorm(6, mean = 4, sd = .3), ncol = 2, nrow = 3, dimnames = list(c('a_i','b_i','c_i')))
colnames(p2) <- c('x','y')
p1 <- matrix(rnorm(6, mean = 2, sd = .6), ncol = 2, nrow = 3, dimnames = list(c('a_i','b_i','c_i')))
colnames(p1) <- c('x','y')
single_point <- matrix(rnorm(1, mean = 2, sd = .6), ncol = 2, nrow = 1, dimnames = list(c('a_i')))
a = c(100,23,300,44)
a1 = c(-90,200,70)
a1_double = c(-90)
a2 = c(20,11,-190,360)
v=c(1); tStep = 0.5

ped_a1 = c(a_n = 20, b_n= 11, C_n = -190, A_n = 360,
       B_n = 33, j_n = -13, k_n = -200, h_n = -45,
       p_n = 320, P_n = -320, l_n = 76)
ped_a2 = c(a_n = 20, b_n= 11, C_n = -190, A_n = 360,
       B_n = 33, j_n = -13, k_n = -200, h_n = -45,
       p_n = 320, P_n = -320, l_n = 76)



testthat::test_that("dist geometry function works", {
  ref = m4ma::dist(p1, p2)
  ll_state = m4ma::dist_rcpp(p1,p2)
  testthat::expect_equal(ll_state, ref)
})

testthat::test_that("dist1 geometry function works", {
  ref = m4ma::dist1(p1 = single_point, p2)
  ll_state = m4ma::dist1_rcpp(p1 = single_point,p2)
  testthat::expect_equal(ll_state, ref)
})

testthat::test_that("angle2 geometry function works", {
  ref = m4ma::angle2(p1, p2)
  ll_state = m4ma::angle2_rcpp(p1,p2)
  testthat::expect_equal(ll_state, ref)
})

testthat::test_that("angle2s geometry function works", {
  ref = m4ma::angle2s(p1,p2)
  ll_state = m4ma::angle2s_rcpp(p1,p2)
  testthat::expect_equal(ll_state, ref)
})

testthat::test_that("aTOd geometry function works", {
  ref = m4ma::aTOd(a)
  ll_state = m4ma::aTOd_rcpp(a)
  testthat::expect_equal(ll_state, ref)
})

testthat::test_that("Iangle geometry function works", {
  ref = m4ma::Iangle(p1,a1,p2)
  ll_state = m4ma::Iangle_rcpp(p1,p2,a1)
  testthat::expect_equal(ll_state, ref)
})

testthat::test_that("Dn geometry function works", {
  ref = m4ma::Dn(p_n = p1, P_n = p2) 
  ll_state = m4ma::Dn_rcpp(p_n = unname(p1), P_n = p2)
  testthat::expect_equal(ll_state, ref)
})

testthat::test_that("minAngle geometry function works", {
  ref = m4ma::minAngle(a1_double, a2) 
  ll_state = m4ma::minAngle_rcpp(a1_double, a2)
  testthat::expect_equal(ll_state, ref)
})

testthat::test_that("scaleVel geometry function works", {
  ref = m4ma::scaleVel(v, tStep) 
  ll_state = m4ma::scaleVel_rcpp(v, tStep)
  testthat::expect_equal(ll_state, ref)
})

testthat::test_that("coneNum geometry function works", {
  ref = m4ma::coneNum(1:33) 
  ll_state = m4ma::coneNum_rcpp(1:33)
  testthat::expect_equal(ll_state, ref)
})

testthat::test_that("ringNum geometry function works", {
  ref = m4ma::ringNum(1:33) 
  ll_state = m4ma::ringNum_rcpp(1:33)
  testthat::expect_equal(ll_state, ref)
})

testthat::test_that("c_vd geometry function works", {
  ref = m4ma::c_vd(cells = 1:33, p1 = c(x=-.25,y=-.25),
                   v1 = v, a1 = a1_double) 
  ll_state = m4ma::c_vd_rcpp(cells = 1:33,p1 = c(x=-.25,y=-.25),
                             v1 = v, a1 = a1_double,
                             vels = matrix(rep(c(1.5, 1, .5), each = 11), ncol = 3),
                             angles = matrix(rep(c(72.5, 50, 32.5, 20, 10, 0, 350, 340, 
                                                   327.5, 310, 287.5), times = 3),ncol = 3))
  testthat::expect_equal(ll_state, ref)
})




