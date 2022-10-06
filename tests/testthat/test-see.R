# Test see functions against reference values

P1 = c(0, 0); P2 = c(0, 1); P3 = c(1, 0); P4 = c(1, 1)


# Line-line intersection
testthat::test_that("Lines intersect", {
  ref = c(0.5, 0.5)
  res = m4ma::line_line_intersection_rcpp(P1, P4, P2, P3)
  testthat::expect_equal(res, ref)
})

testthat::test_that("Lines do not intersect", {
  ref = c(Inf, Inf)
  res = m4ma::line_line_intersection_rcpp(P1, P2, P3, P4)
  testthat::expect_equal(res, ref)
})

# If two lines intersect outside of their spans
testthat::test_that("Lines intersect outside", {
  ref = as.numeric(c(NA, NA))
  res = m4ma::line_line_intersection_rcpp(
    P1, c(0, 2), P2, P4, interior_only = TRUE
  )
  testthat::expect_equal(res, ref)
})


# Sees goal
testthat::test_that("Does not see goal", {
  objects = list(
    list(x = c(-0.25, 0.75), y = c(0.25, 0.75)),
    list(x = c(0.75, 0.75), y = c(0.25, 0.25))
  )
  res = m4ma::seesGoal_rcpp(P1, P4, objects)
  testthat::expect_false(res)
})

testthat::test_that("Sees goal", {
  objects = list(
    list(x = c(0.25, 0.75), y = c(0.25, 0.75)),
    list(x = c(0.75, 0.75), y = c(0.25, 0.25))
  )
  res = m4ma::seesGoal_rcpp(P1, P4, objects)
  testthat::expect_true(res)
})

state = list(
  p = matrix(P1, 1, 2),
  P = list(
    matrix(c(0.5, 0.5), 1, 2)
  )
)


# Sees current goal
testthat::test_that("Does not see current goal", {
  state = list(
    p = matrix(P1, 1, 2),
    P = list(
      matrix(P4, 1, 2)
    )
  )
  attr(state$P[[1]], "i") = 1
  objects = list(
    list(x = c(-0.25, 0.75), y = c(0.25, 0.75)),
    list(x = c(0.75, 0.75), y = c(0.25, 0.25))
  )
  res = m4ma::seesCurrentGoal_rcpp(1, state, objects)
  testthat::expect_false(res)
})

testthat::test_that("Sees current goal", {
  state = list(
    p = matrix(P1, 1, 2),
    P = list(
      matrix(P4, 1, 2)
    )
  )
  attr(state$P[[1]], "i") = 1
  objects = list(
    list(x = c(0.25, 0.75), y = c(0.25, 0.75)),
    list(x = c(0.75, 0.75), y = c(0.25, 0.25))
  )
  res = m4ma::seesCurrentGoal_rcpp(1, state, objects)
  testthat::expect_true(res)
})


# Sees many goals
testthat::test_that("Does not see any goal", {
  ref = logical(3)
  ps = rbind(P4, c(2, 2), c(3, 3))
  objects = list(
    list(x = c(-0.25, 0.75), y = c(0.25, 0.75)),
    list(x = c(0.75, 0.75), y = c(0.25, 0.25))
  )
  res = m4ma::seesMany_rcpp(P1, ps, objects)
  testthat::expect_equal(res, ref)
})

testthat::test_that("Sees some goals", {
  ref = c(FALSE, TRUE, FALSE)
  ps = rbind(P2, P3, P4)
  objects = list(
    list(x = c(-0.25, 0.75), y = c(0.25, 0.75)),
    list(x = c(0.75, 0.75), y = c(0.25, 0.25))
  )
  res = m4ma::seesMany_rcpp(P1, ps, objects)
  testthat::expect_equal(res, ref)
})

testthat::test_that("Sees all goals", {
  ref = c(TRUE, TRUE, TRUE)
  ps = rbind(c(2, 0), P3, c(3, 0))
  objects = list(
    list(x = c(-0.25, 0.75), y = c(0.25, 0.75)),
    list(x = c(0.75, 0.75), y = c(0.25, 0.25))
  )
  res = m4ma::seesMany_rcpp(P1, ps, objects)
  testthat::expect_equal(res, ref)
})


# Sees goals from ok cells
testthat::test_that("Sees goal from ok cells", {
  ref = as.logical(c(
    rep(0, 5),
    rep(1, 3),
    0, 1, 0, 1, 1,
    rep(0, 3),
    1,
    rep(0, 4),
    rep(1, 6),
    rep(0, 3),
    1, 0, 0
  ))
  state = list(
    p = matrix(P1, 1, 2),
    P = list(
      matrix(P4, 1, 2)
    )
  )
  attr(state$P[[1]], "i") = 1
  objects = list(
    list(x = c(-0.25, 0.75), y = c(0.25, 0.75)),
    list(x = c(0.75, 0.75), y = c(0.25, 0.25))
  )
  set.seed(123123)
  centres = matrix(rnorm(66), 33, 2)
  ok = as.logical(c(0, rep(1, 32)))
  res = m4ma::seesGoalOK_rcpp(1, objects, state, centres, ok)
  testthat::expect_equal(res, ref)
})
