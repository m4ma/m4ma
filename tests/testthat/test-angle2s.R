set.seed(123)
p2 <- matrix(rnorm(6, mean = 4, sd = .3), ncol = 2, nrow = 3, dimnames = list(c('a_i','b_i','c_i')))
colnames(p2) <- c('x','y')
p1 <- matrix(rnorm(6, mean = 2, sd = .6), ncol = 2, nrow = 3, dimnames = list(c('a_i','b_i','c_i')))
colnames(p1) <- c('x','y')

testthat::test_that("angle2s works", {
  ref = m4ma::angle2s(p1,p2)
  ll_state = m4ma::angle2s_rcpp(p1,p2)
  testthat::expect_equal(ll_state, ref)
})

