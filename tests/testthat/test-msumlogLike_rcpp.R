
load('tests/data/trace_i.rda')

nests = list(
  Central = c(0, 6, 17, 28),
  NonCentral = c(0:33)[-c(6, 17, 28)],
  acc = c(1:11),
  const = c(12:22),
  dec = c(0, 23:33)
)

alpha = list(
  Central = rep(1/3, 4),
  NonCentral = c(1/3, rep(0.5, 4), 1/3, rep(0.5, 9), 1/3, 
                 rep(0.5, 9), 1/3, rep(0.5, 5)),
  acc = c(rep(0.5, 4), 1, 1/3, rep(0.5, 5)),
  const = c(rep(0.5, 4), 1, 1/3, rep(0.5, 5)),
  dec = c(1/3, rep(0.5, 4), 1, 1/3, rep(0.5, 5))
)

cell_nest = cbind(
  t(matrix(c(c(1, 5), rep(2:3, 5), c(1, 3), rep(2:3, 5), rep(c(2, 4), 5), 
             c(1, 4), rep(c(2, 4), 5), rep(c(2, 5), 5), c(1, 5), 
             rep(c(2, 5), 5)),
           nrow = 2)),
  cbind(c(1, 1:5, 2, 6:15, 3, 16:25, 4, 26:30), 
        c(1, 1:11, 1:11, 2:12)))-1

p = attr(trace_i, 'pMat')

testthat::test_that("State likelihood works", {
  ref = 0.5
  ll_state = m4ma::like_state_rcpp(trace_i[[1]], p, n = 3, nests, alpha, m4ma::get_cell_nest())
  testthat::expect_equal(ll_state, ref)
})

testthat::test_that("State likelihood loop works", {
  ref = c(-Inf, -487.1970260, -0.6931472, -Inf, -14.8647133, -Inf, -Inf, -309.3628411, -Inf)
  ll_states = m4ma::like_states_rcpp(trace_i, p, nests, alpha, m4ma::get_cell_nest())
  testthat::expect_equal(log(unlist(ll_states)), ref)
})

testthat::test_that("Trace likelihood sum works", {
  ref = 176.738817
  ll_sum = m4ma::msumlogLike_rcpp(p, trace_i, minLike = 1e-10, mult = -1)
  testthat::expect_equal(ll_sum, ref)
})
