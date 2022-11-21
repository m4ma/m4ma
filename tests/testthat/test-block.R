# Test rcpp block functions against R counterparts

# Create dummy objects for testing

objects = list(
  list(x = c(0.15, 0.35), y = c(0.35, 0.55)),
  list(x = c(0.55, 0.35), y = c(0.75, 0.65))
)

set.seed(123)

centres = matrix(rnorm(66), 33, 2)

r = 0.5

ok = as.logical(rbinom(33, size = 1, prob = 0.5))

testthat::test_that("object2lines function works", {
  ref = m4ma::object2lines_r(objects[[1]])
  res = m4ma::object2lines_rcpp(objects[[1]])
  testthat::expect_equal(res, list(P1 = ref[,,1], P2 = ref[,,2]))
})

testthat::test_that("bodyObjectOverlap function works", {
  ok_centres = centres[ok, , drop = FALSE]
  ref = m4ma::bodyObjectOverlap_r(object2lines_r(objects[[1]]), r, ok_centres)
  res = m4ma::bodyObjectOverlap_rcpp(object2lines_rcpp(objects[[1]]), r,
                                     ok_centres)
  testthat::expect_equal(res, ref)
})

testthat::test_that("bodyObjectOK function works", {
  ref = m4ma::bodyObjectOK_r(r, centres, objects, ok)
  res = m4ma::bodyObjectOK_rcpp(r, centres, objects, ok)
  testthat::expect_equal(res, ref)
})
