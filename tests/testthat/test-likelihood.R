# Test likelihood functions against reference values on test trace

test_filepath = file.path('data', 'trace_i.rda')

test_obj_name = load(test_filepath)

p = attr(get(test_obj_name), 'pMat')

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

test_trace_rcpp = m4ma::create_rcpp_trace(get(test_obj_name))

ref = -176.738817

testthat::test_that("Subject likelihood computation works", {
  ll_subject = sapply(1:length(test_trace_rcpp), function(i) {
    sum(sapply(1:length(test_trace_rcpp[[i]]), function(j) {
      m4ma::like_subject(test_trace_rcpp[[i]][[j]], p[j, ], n = j - 1, nests, alpha, m4ma::get_cell_nest())
    }))
  })
  testthat::expect_equal(sum(ll_subject), ref)
})

testthat::test_that("State likelihood computation works", {
  ll_state = sapply(test_trace_rcpp, function(state) {
    m4ma::like_state(state, p, nests, alpha, m4ma::get_cell_nest())
  })
  testthat::expect_equal(sum(ll_state), ref)
})

testthat::test_that("Trace likelihood sum works", {
  ll_sum = m4ma::msumlogLike(p, test_trace_rcpp, nests, alpha, m4ma::get_cell_nest())
  testthat::expect_equal(-ll_sum, ref)
})
