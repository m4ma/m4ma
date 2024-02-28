# Test rcpp utility functions against R counterparts on test trace

test_filepath = testthat::test_path('data', 'trace_i.rda')

test_obj_name = load(test_filepath)

p = exp(attr(get(test_obj_name), 'pMat'))


testthat::test_that("Blocked-angle utility computation works (all cells)", {
  BA <- runif(33, 0, pi)
  
  names(BA) <- as.character(1:33)
  
  ref = m4ma::baUtility_r(c(bBA = 4, aBA = 2), BA)
  ba_utility = m4ma::baUtility_rcpp(2, 4, BA, 0:32)
  testthat::expect_equal(ba_utility, ref)
})

testthat::test_that("Blocked-angle utility computation works (subset cells)", {
  BA <- runif(10, 0, pi)
  BA_idx <- sample(1:33, 10)
  
  names(BA) <- as.character(BA_idx)
  
  ref = m4ma::baUtility_r(c(bBA = 4, aBA = 2), BA)
  ba_utility = m4ma::baUtility_rcpp(2, 4, BA, BA_idx-1)
  testthat::expect_equal(ba_utility, ref)
})

testthat::test_that("Blocked-angle utility computation works (exponential scale)", {
  BA <- runif(33, 0, pi)
  
  names(BA) <- as.character(1:33)
  
  ref = m4ma::baUtility_r(exp(c(bBA = 4, aBA = 2)), BA)
  ba_utility = m4ma::baUtility_rcpp(exp(2), exp(4), BA, 0:32)
  testthat::expect_equal(ba_utility, ref)
})

testthat::test_that("Current-angle utility computation works", {
  ref = unname(m4ma::caUtility_r(c(aCA = 1, bCA = 2, bCAlr = 10)))
  ca_utility = m4ma::caUtility_rcpp(1, 2, 10)
  testthat::expect_equal(ca_utility, ref)
})

testthat::test_that("Current-angle utility computation works (exponential scale)", {
  ref = unname(m4ma::caUtility_r(exp(c(aCA = 1, bCA = 2, bCAlr = 10))))
  ca_utility = m4ma::caUtility_rcpp(exp(1), exp(2), exp(10))
  testthat::expect_equal(ca_utility, ref)
})

testthat::test_that("Follow-leader utility computation works", {
  leaders = matrix(c(0, 0.5, 0, 6, 1, 0), 3, 2, dimnames = list(c("cell", "angleDisagree", "inGroup"), NULL))
  dists = matrix(runif(66), 2, 33)
  
  ref = unname(m4ma::flUtility_r(c(aFL = 2, bFL = 1, dFL = 0), list(leaders = leaders, dists = dists)))
  fl_utility = m4ma::flUtility_rcpp(2, 1, 0, leaders, dists)
  testthat::expect_equal(fl_utility, ref)
})


testthat::test_that("Follow-leader utility computation works (exponential scale)", {
  leaders = matrix(c(0, 0.5, 0, 6, 1, 0), 3, 2, dimnames = list(c("cell", "angleDisagree", "inGroup"), NULL))
  dists = matrix(runif(66), 2, 33)
  
  ref = unname(m4ma::flUtility_r(exp(c(aFL = 2, bFL = 1, dFL = 0)), list(leaders = leaders, dists = dists)))
  fl_utility = m4ma::flUtility_rcpp(exp(2), exp(1), exp(0), leaders, dists)
  testthat::expect_equal(fl_utility, ref)
})


testthat::test_that("Goal-angle utility computation works", {
  GA <- runif(11, 0, pi)
  
  ref = m4ma::gaUtility_r(c(bGA = 10, aGA = 2), GA)
  ga_utility = m4ma::gaUtility_rcpp(10, 2, GA)
  testthat::expect_equal(ga_utility, ref)
})


testthat::test_that("Goal-angle utility computation works", {
  GA <- runif(11, 0, pi)
  
  ref = m4ma::gaUtility_r(exp(c(bGA = 10, aGA = 2)), GA)
  ga_utility = m4ma::gaUtility_rcpp(exp(10), exp(2), GA)
  testthat::expect_equal(ga_utility, ref)
})


testthat::test_that("Interpersonal distance utility computation works", {
  ID <- matrix(runif(66, 0, 5), 2, 33, dimnames = list(c("a1", "a2"), NULL))
  ok <- matrix(as.logical(rbinom(33, 1, 0.2), 11, 3))
  group <- c(a0 = 1, a1 = 1, a2 = 2)
  n <- 1
  u <- as.vector(ifelse(ok, 0, -Inf))
  namesInGroup <- names(group[-n][group[-n] == group[n]])
  inGroup <- dimnames(ID)[[1]] %in% namesInGroup
  
  ref = m4ma::idUtility_r(c(bID = 1, aID = 2, dID = 0), n, ID, ok, group)
  id_utility = m4ma::idUtility_rcpp(1, 0, 2, inGroup, ok, ID, u)
  testthat::expect_equal(id_utility, ref)
})


testthat::test_that("Interpersonal distance utility computation works (exponential scale)", {
  ID <- matrix(runif(66, 0, 5), 2, 33, dimnames = list(c("a1", "a2"), NULL))
  ok <- matrix(as.logical(rbinom(33, 1, 0.2), 11, 3))
  group <- c(a0 = 1, a1 = 1, a2 = 2)
  n <- 1
  u <- as.vector(ifelse(ok, 0, -Inf))
  namesInGroup <- names(group[-n][group[-n] == group[n]])
  inGroup <- dimnames(ID)[[1]] %in% namesInGroup
  
  ref = m4ma::idUtility_r(exp(c(bID = 1, aID = 2, dID = 0)), n, ID, ok, group)
  id_utility = m4ma::idUtility_rcpp(exp(1), exp(0), exp(2), inGroup, ok, ID, u)
  testthat::expect_equal(id_utility, ref)
})


testthat::test_that("Interpersonal distance utility computation works (bID == 0)", {
  ID <- matrix(runif(66, 0, 5), 2, 33, dimnames = list(c("a1", "a2"), NULL))
  ok <- matrix(as.logical(rbinom(33, 1, 0.2), 11, 3))
  group <- c(a0 = 1, a1 = 1, a2 = 2)
  n <- 1
  u <- as.vector(ifelse(ok, 0, -Inf))
  namesInGroup <- names(group[-n][group[-n] == group[n]])
  inGroup <- dimnames(ID)[[1]] %in% namesInGroup
  
  ref = m4ma::idUtility_r(c(bID = 0, aID = 2, dID = 0), n, ID, ok, group)
  id_utility = m4ma::idUtility_rcpp(0, 0, 2, inGroup, ok, ID, u)
  testthat::expect_equal(id_utility, ref)
})


testthat::test_that("Interpersonal distance utility computation works (all not ok)", {
  ID <- matrix(runif(66, 0, 5), 2, 33, dimnames = list(c("a1", "a2"), NULL))
  ok <- matrix(as.logical(FALSE, 11, 3))
  group <- c(a0 = 1, a1 = 1, a2 = 2)
  n <- 1
  u <- as.vector(ifelse(ok, 0, -Inf))
  namesInGroup <- names(group[-n][group[-n] == group[n]])
  inGroup <- dimnames(ID)[[1]] %in% namesInGroup
  
  ref = m4ma::idUtility_r(c(bID = 1, aID = 2, dID = 0), n, ID, ok, group)
  id_utility = m4ma::idUtility_rcpp(1, 0, 2, inGroup, ok, ID, u)
  testthat::expect_equal(id_utility, ref)
})


testthat::test_that("Preferred speed utility computation works", {
  v <- 1
  d <- 5
  
  ref = unname(m4ma::psUtility_r(c(aPS = 2, bPS = 4, sPref = 1, sSlow = 1), v, d))
  ps_utility = m4ma::psUtility_rcpp(2, 4, 1, 1, v, d)
  testthat::expect_equal(ps_utility, ref)
})


testthat::test_that("Preferred speed utility computation works", {
  v <- 1
  d <- 5
  
  ref = unname(m4ma::psUtility_r(exp(c(aPS = 2, bPS = 4, sPref = 1, sSlow = 1)), v, d))
  ps_utility = m4ma::psUtility_rcpp(exp(2), exp(4), exp(1), exp(1), v, d)
  testthat::expect_equal(ps_utility, ref)
})


testthat::test_that("Walk-beside utility computation works", {
  buddies = matrix(c(0, 0.5, 6, 1), 2, 2, dimnames = list(c("cell", "angleDisagree"), NULL))
  dists = matrix(runif(66), 2, 33)
  
  ref = unname(m4ma::wbUtility_r(c(bWB = 2, aWB = 2),list(buddies = buddies, dists = dists)))
  wb_utility = m4ma::wbUtility_rcpp(2, 2, buddies, dists)
  testthat::expect_equal(wb_utility, ref)
})


testthat::test_that("Walk-beside utility computation works (exponential scale)", {
  buddies = matrix(c(0, 0.5, 6, 1), 2, 2, dimnames = list(c("cell", "angleDisagree"), NULL))
  dists = matrix(runif(66), 2, 33)
  
  ref = unname(m4ma::wbUtility_r(exp(c(bWB = 2, aWB = 2)),list(buddies = buddies, dists = dists)))
  wb_utility = m4ma::wbUtility_rcpp(exp(2), exp(2), buddies, dists)
  testthat::expect_equal(wb_utility, ref)
})


testthat::test_that("Utility computation works", {
  ref = c(-19915374, -29604211, -31746385)

  u = sapply(1:length(get(test_obj_name)), function(j) {
    state = get(test_obj_name)[[j]]
    out = sapply(1:length(state$v), function(i) {
      group <- state$group
      ok <- state$ok[[i]]
      if (!is.null(state$ID[[i]])) {
        namesInGroup <- names(group[-i][group[-i] == group[i]])
        inGroup <- dimnames(state$ID[[i]])[[1]] %in% namesInGroup
        ok <- ok & apply(state$ID[[i]], 2, function(x) {
          all(x > 0)
        })
      } else {
        inGroup <- FALSE
      }
      
      u = m4ma::utility(
        p[i, ], state$v[i], state$d[i], pmax(state$BA[[i]], 0), state$GA[[i]],
        state$ID[[i]], inGroup, as.vector(ifelse(ok, 0, -Inf)), state$FL[[i]], state$WB[[i]], ok
      )

      return(sum(u[is.finite(u)]))
    })

    return(sum(out))
  })

  testthat::expect_equal(u, ref)
})
