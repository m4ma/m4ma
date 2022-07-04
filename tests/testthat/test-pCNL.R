
n_cells = 34

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

mu = 1
mum = rep(1, 5)

testthat::test_that("Probability computation works", {
  ref = rep(0, n_cells)
  ref[12:13] = 0.5
  set.seed(123556)
  utility = runif(n_cells, -10e04, 0)
  prob = sapply(1:n_cells, function(i) m4ma::pcnl_rcpp(cell_nest[i, ], utility, mum, nests, alpha, mu))
  testthat::expect_equal(prob, ref)
})

testthat::test_that("Infinite utility handling works", {
  ref = rep(0, n_cells)
  ref[1:2] = c(2, 1)/3
  set.seed(123556)
  utility = runif(n_cells, -10e04, 0)
  utility[1] = Inf
  prob = sapply(1:n_cells, function(i) m4ma::pcnl_rcpp(cell_nest[i, ], utility, mum, nests, alpha, mu))
  testthat::expect_equal(prob, ref)
})

testthat::test_that("All infinite utility handling works", {
  utility_inf = rep(Inf, n_cells)
  prob_inf = sapply(1:n_cells, function(i) m4ma::pcnl_rcpp(cell_nest[i, ], utility_inf, mum, nests, alpha, mu))
  utility_minf = rep(-Inf, n_cells)
  prob_minf = sapply(1:n_cells, function(i) m4ma::pcnl_rcpp(cell_nest[i, ], utility_minf, mum, nests, alpha, mu))
  testthat::expect_equal(prob_inf, prob_minf)
})

testthat::test_that("Infinite alpha handling works", {
  ref = c(1, rep(0, n_cells-1))
  alpha_inf = alpha
  alpha_inf[["Central"]][1] = Inf
  set.seed(123556)
  utility = runif(n_cells, -10e04, 0)
  prob = sapply(1:n_cells, function(i) m4ma::pcnl_rcpp(cell_nest[i, ], utility, mum, nests, alpha_inf, mu))
  testthat::expect_equal(prob, ref)
})

testthat::test_that("Infinite mu handling works", {
  ref = rep(0, n_cells)
  mu_inf = Inf
  set.seed(123556)
  utility = runif(n_cells, -10e04, 0)
  prob = sapply(1:n_cells, function(i) m4ma::pcnl_rcpp(cell_nest[i, ], utility, mum, nests, alpha, mu_inf))
  testthat::expect_equal(prob, ref)
})

testthat::test_that("Infinite mum handling works", {
  ref = rep(0, n_cells)
  ref[12:13] = c(0.25, 0.25)
  mum_inf = c(Inf, rep(1, 4))
  set.seed(123556)
  utility = runif(n_cells, -10e04, 0)
  prob = sapply(1:n_cells, function(i) m4ma::pcnl_rcpp(cell_nest[i, ], utility, mum_inf, nests, alpha, mu))
  testthat::expect_equal(prob, ref)
})
