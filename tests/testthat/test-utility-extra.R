# Test utility helper functions against R counterparts

# Pedestrian names
nms = c('a1', 'a2', 'c1')

# Current pedestrian
p1 = matrix(c(0, 0), 1, 2)
rownames(p1) = nms[1]

# All pedestrians
p2 = rbind(c(0, 0), c(0.5, 0.5), c(1, 1))
rownames(p2) = nms

# Pedestrian angle in front
a_front = 90
names(a_front) = nms[1]

# Pedestrian angles not in front
a_not_front = 135
names(a_not_front) = nms[1]

# Pedestrian angles vector
a = c(90, 90, 90)
names(a) = nms

# Pedestrian velocities
v = c(0.5, 0.5, 0.5)
names(v) = nms

# Goals
P1 = matrix(c(2, 2, -1, -1), 2, 2)
rownames(P1) = c('g1', 'g2')
attr(P1, 'i') = 1

# Pedestrian radius
r = c(0.5, 0.5, 0.5)
names(r) = nms

# Pedestrian group
group = c(2, 1, 2)
names(group) = nms

# Cell centres
set.seed(23123)
centres = matrix(rnorm(66), 33, 2)

# Predicted pedestrian positions
p_pred = rbind(c(0, 0), c(0.75, 0.75), c(1.25, 1.25))
rownames(p_pred) = nms

# Objects not occluding view between p1 and p2
objects = list(
  list(x = c(1, 0), y = c(0, 1)),
  list(x = c(0, 1), y = c(1, 0))
)

# Objects occluding view between p1 and p2
objects_occlude = list(
  list(x = c(0.15, 0.35), y = c(0.15, 0.35))
)

# State list
state = list(
  p = p2,
  a = c(45, 45, 45),
  v = v,
  r = r,
  P = list(
    P1, P1, P1
  ),
  group = group
)
names(state$a) = nms

# Pos of current pedestrian in p2
n = 1

testthat::test_that("Destination angle computation works", {
  ref = m4ma::destinationAngle_r(a_front, p1, P1)
  res = m4ma::destinationAngle_rcpp(a_front, p1, P1)
  testthat::expect_equal(res, ref)
})

testthat::test_that(
  "Predicted distance works (when all pedestrians are seen in front)", {
    ref = m4ma::predClose_r(n, p1, a_front,
                            p2, r, centres, p_pred, objects)
    res = m4ma::predClose_rcpp(n, p1, a_front,
                               p2, r, centres, p_pred, objects)
    testthat::expect_equal(res, ref)
  }
)

testthat::test_that(
  "Predicted distance works (when one or less pedestrians are occluded)", {
    ref = m4ma::predClose_r(n, p1, a_front,
                            p2, r, centres, p_pred, objects_occlude)
    res = m4ma::predClose_rcpp(n, p1, a_front,
                               p2, r, centres, p_pred, objects_occlude)
    testthat::expect_equal(res, ref)
  }
)

testthat::test_that(
  "Predicted distance is NULL (when one or less pedestrians are occluded)", {
    res = m4ma::predClose_rcpp(n, p1, a_front,
                               p2, r, centres, p_pred, objects_occlude)
    testthat::expect_equal(res, NULL)
  }
)

testthat::test_that(
  "Predicted distance works (when pedestrians are seen but not in front)", {
    ref = m4ma::predClose_r(n, p1, a_not_front,
                            p2, r, centres, p_pred, objects)
    res = m4ma::predClose_rcpp(n, p1, a_not_front,
                               p2, r, centres, p_pred, objects_occlude)
    testthat::expect_equal(res, ref)
  }
)

testthat::test_that(
  "Predicted distance is NULL (when pedestrians are seen but not in front)", {
    res = m4ma::predClose_rcpp(n, p1, a_not_front,
                               p2, r, centres, p_pred, objects_occlude)
    testthat::expect_equal(res, NULL)
  }
)

testthat::test_that("Egocentric objects computation works", {
  ref = m4ma::eObjects_r(p1, p2, r)
  ref_list = list(ac = ref[, , 1], cw = ref[, , 2]) # C++ requires list instead of array
  res = m4ma::eObjects_rcpp(p1, p2, r)
  testthat::expect_equal(res, ref_list)
})

testthat::test_that(
  "Intersecting cones works (when all pedestrians are seen in front)", {
  ref = m4ma::iCones_r(p1, a_front, p2[-n, , drop = FALSE], r, objects)
  res = m4ma::iCones_rcpp(p1, a_front, p2[-n, , drop = FALSE], r, objects)
  testthat::expect_equal(res, ref)
})

testthat::test_that(
  "Intersecting cones works (when single pedestrian)", {
    ref = m4ma::iCones_r(p1, a_front, matrix(, 0, 2), r, objects)
    res = m4ma::iCones_rcpp(p1, a_front, matrix(, 0, 2), r, objects)
    testthat::expect_equal(res, ref)
  }
)

testthat::test_that(
  "Intersecting cones is NULL (when single pedestrian)", {
    res = m4ma::iCones_rcpp(p1, a_front, matrix(, 0, 2), r, objects)
    testthat::expect_equal(res, NULL)
  }
)

testthat::test_that(
  "Intersecting cones works (when one other pedestrian)", {
    ref = m4ma::iCones_r(p1, a_front, p2[2, , drop = FALSE], r, objects)
    res = m4ma::iCones_rcpp(p1, a_front, p2[2, , drop = FALSE], r, objects)
    testthat::expect_equal(res, ref)
  }
)

testthat::test_that(
  "Intersecting cones works (when no intersection)", {
    ref = m4ma::iCones_r(p1, 175, p2[-n, , drop = FALSE], r, objects)
    res = m4ma::iCones_rcpp(p1, 175, p2[-n, , drop = FALSE], r, objects)
    testthat::expect_equal(res, ref)
  }
)

testthat::test_that(
  "Intersecting cones is NULL (when no intersection)", {
    res = m4ma::iCones_rcpp(p1, 175, p2[-n, , drop = FALSE], r, objects)
    testthat::expect_equal(res, NULL)
  }
)

testthat::test_that(
  "Intersecting cones works (when one cone intersecting and seen)", {
    ref = m4ma::iCones_r(p1, a_not_front, p2[3, , drop = FALSE], r, objects)
    res = m4ma::iCones_rcpp(p1, a_not_front, p2[3, , drop = FALSE], r, objects)
    testthat::expect_equal(res, ref)
  }
)

testthat::test_that(
  "Intersecting cones works (when one cone intersecting but not seen)", {
    ref = m4ma::iCones_r(p1, a_not_front,
                         p2[3, , drop = FALSE], r, objects_occlude)
    res = m4ma::iCones_rcpp(p1, a_not_front,
                            p2[3, , drop = FALSE], r, objects_occlude)
    testthat::expect_equal(res, ref)
  }
)

testthat::test_that(
  "Intersecting cones is NULL (when one cone intersecting but not seen)", {
    res = m4ma::iCones_rcpp(p1, a_not_front,
                            p2[3, , drop = FALSE], r, objects_occlude)
    testthat::expect_equal(res, NULL)
  }
)

testthat::test_that(
  "Intersecting cones works (when all pedestrians are occluded in front)", {
    ref = m4ma::iCones_r(p1, a_front,
                         p2[-n, , drop = FALSE], r, objects_occlude)
    res = m4ma::iCones_rcpp(p1, a_front,
                            p2[-n, , drop = FALSE], r, objects_occlude)
    testthat::expect_equal(res, ref)
  })

testthat::test_that("Intersecting cones to cells transformation works", {
  # Make test independent of iCones
  cones = c(rnorm(6))
  names(cones) = as.character(1:6)
  v = 0.5
  ref = m4ma::iCones2Cells_r(cones, v)
  res = m4ma::iCones2Cells_rcpp(cones, v)
  testthat::expect_equal(res, ref)
})

testthat::test_that(
  "Blocked angle works (when intersecting cones and seen)", {
    ref = m4ma::blockedAngle_r(n, state, p_pred, objects)
    res = m4ma::blockedAngle_rcpp(n, state, p_pred, objects)
    testthat::expect_equal(res, if (length(ref) == 0) unname(ref) else ref)
  }
)

testthat::test_that(
  "Blocked angle works (when no intersecting cones)", {
    a_back = c(175, 175, 175)
    names(a_back) = nms
    state_occluded = state
    state_occluded$a = a_back
    ref = m4ma::blockedAngle_r(n, state_occluded, p2, objects_occlude)
    res = m4ma::blockedAngle_rcpp(n, state_occluded, p_pred, objects_occlude)
    testthat::expect_equal(res, if (length(ref) == 0) unname(ref) else ref)
  }
)

testthat::test_that(
  "Blocked angle returns numeric(0) (when no intersecting cones)", {
    a_back = c(175, 175, 175)
    names(a_back) = nms
    state_occluded = state
    state_occluded$a = a_back
    res = m4ma::blockedAngle_rcpp(n, state_occluded, p_pred, objects_occlude)
    testthat::expect_equal(res, numeric(0))
  }
)

testthat::test_that(
  "Get leaders works (when seen and in front)", {
    state_leaders = state
    ref = m4ma::getLeaders_r(n, state_leaders, centres, objects)
    res = m4ma::getLeaders_rcpp(n, state_leaders, centres, objects)
    testthat::expect_equal(res, ref)
  }
)

testthat::test_that(
  "Get leaders works (when not seen and in front)", {
    state_leaders = state
    ref = m4ma::getLeaders_r(n, state_leaders, centres, objects_occlude)
    res = m4ma::getLeaders_rcpp(n, state_leaders, centres, objects_occlude)
    testthat::expect_equal(res, ref)
  }
)

testthat::test_that(
  "Get leaders works (when seen and not in front)", {
    state_leaders = state
    state_leaders$a = a
    ref = m4ma::getLeaders_r(n, state_leaders, centres, objects)
    res = m4ma::getLeaders_rcpp(n, state_leaders, centres, objects)
    testthat::expect_equal(res, ref)
  }
)

testthat::test_that(
  "Get leaders returns NULL (when seen and not in front)", {
    state_leaders = state
    state_leaders$a = a
    res = m4ma::getLeaders_rcpp(n, state_leaders, centres, objects)
    testthat::expect_equal(res, NULL)
  }
)

testthat::test_that(
  "Get leaders works (when seen and not in cones)", {
    state_leaders = state
    p2_back = rbind(c(0, 0), c(-0.5, -0.5), c(-1, -1))
    rownames(p2_back) = nms
    state_leaders$p = p2_back
    ref = m4ma::getLeaders_r(n, state_leaders, centres, objects)
    res = m4ma::getLeaders_rcpp(n, state_leaders, centres, objects)
    testthat::expect_equal(res, ref)
  }
)

testthat::test_that(
  "Get leaders returns NULL (when seen and not in cones)", {
    state_leaders = state
    p2_back = rbind(c(0, 0), c(-0.5, -0.5), c(-1, -1))
    rownames(p2_back) = nms
    state_leaders$p = p2_back
    res = m4ma::getLeaders_rcpp(n, state_leaders, centres, objects)
    testthat::expect_equal(res, NULL)
  }
)

testthat::test_that(
  "Get leaders works (when no in-group and only group)", {
    state_leaders = state
    group = c(1, 2, 3)
    names(group) = nms
    state_leaders$group = group
    ref = m4ma::getLeaders_r(n, state_leaders, centres,
                             objects, onlyGroup = TRUE)
    res = m4ma::getLeaders_rcpp(n, state_leaders, centres,
                                objects, onlyGroup = TRUE)
    testthat::expect_equal(res, ref)
  }
)

testthat::test_that(
  "Get leaders returns NULL (when no in-group and only group)", {
    state_leaders = state
    group = c(1, 2, 3)
    names(group) = nms
    state_leaders$group = group
    res = m4ma::getLeaders_rcpp(n, state_leaders, centres,
                                objects, onlyGroup = TRUE)
    testthat::expect_equal(res, NULL)
  }
)

testthat::test_that(
  "Get leaders works (when no in-group)", {
    state_leaders = state
    group = c(1, 2, 3)
    names(group) = nms
    state_leaders$group = group
    ref = m4ma::getLeaders_r(n, state_leaders, centres,
                             objects)
    res = m4ma::getLeaders_rcpp(n, state_leaders, centres,
                                objects)
    testthat::expect_equal(res, ref)
  }
)

testthat::test_that(
  "Get leaders works (when no in-group and no prefer group)", {
    state_leaders = state
    group = c(1, 2, 3)
    names(group) = nms
    state_leaders$group = group
    ref = m4ma::getLeaders_r(n, state_leaders, centres,
                             objects, preferGroup = FALSE)
    res = m4ma::getLeaders_rcpp(n, state_leaders, centres,
                                objects, preferGroup = FALSE)
    testthat::expect_equal(res, ref)
  }
)

testthat::test_that(
  "Get leaders works (when different directions)", {
    state_leaders = state
    a = c(0, 175, 175)
    names(a) = nms
    state_leaders$a = a
    ref = m4ma::getLeaders_r(n, state_leaders, centres, objects)
    res = m4ma::getLeaders_rcpp(n, state_leaders, centres, objects)
    testthat::expect_equal(res, ref)
  }
)

testthat::test_that(
  "Get leaders returns NULL (when different directions)", {
    state_leaders = state
    a = c(0, 175, 175)
    names(a) = nms
    state_leaders$a = a
    res = m4ma::getLeaders_rcpp(n, state_leaders, centres, objects)
    testthat::expect_equal(res, NULL)
  }
)

testthat::test_that(
  "Get leaders works (when duplicate leaders)", {
    state_leaders = state
    group = c(2, 2, 2)
    names(group) = nms
    state_leaders$group = group
    p2_same = rbind(c(0, 0), c(0.5, 0), c(0, 0.5))
    rownames(p2_same) = c(nms[1], nms[2], nms[2])
    state_leaders$p = p2_same
    ref = m4ma::getLeaders_r(n, state_leaders, centres, objects)
    res = m4ma::getLeaders_rcpp(n, state_leaders, centres, objects)
    testthat::expect_equal(res, ref)
  }
)

testthat::test_that(
  "Get leaders works (when pick best)", {
    state_leaders = state
    group = c(2, 2, 2)
    names(group) = nms
    state_leaders$group = group
    p2_same = rbind(c(0, 0), c(0.5, 0), c(0, 0.5))
    rownames(p2_same) = nms
    state_leaders$p = p2_same
    ref = m4ma::getLeaders_r(n, state_leaders, centres,
                             objects, pickBest = TRUE)
    res = m4ma::getLeaders_rcpp(n, state_leaders, centres,
                                objects, pickBest = TRUE)
    testthat::expect_equal(res, ref)
  }
)

testthat::test_that(
  "Get buddy works (when seen and in front)", {
    state_buddy = state
    ref = m4ma::getBuddy_r(n, group, a, p_pred, centres,
                           objects, FALSE, state_buddy)
    res = m4ma::getBuddy_rcpp(n, group, a, p_pred, centres,
                           objects, FALSE, state_buddy)
    testthat::expect_equal(res, ref)
  }
)

testthat::test_that(
  "Get buddy works (when not seen and in front)", {
    state_buddy = state
    ref = m4ma::getBuddy_r(n, group, a, p_pred, centres,
                           objects_occlude, FALSE, state_buddy)
    res = m4ma::getBuddy_rcpp(n, group, a, p_pred, centres,
                           objects_occlude, FALSE, state_buddy)
    testthat::expect_equal(res, ref)
  }
)

testthat::test_that(
  "Get buddy returns NULL (when not seen and in front)", {
    state_buddy = state
    res = m4ma::getBuddy_rcpp(n, group, a, p_pred, centres,
                              objects_occlude, FALSE, state_buddy)
    testthat::expect_equal(res, NULL)
  }
)

testthat::test_that(
  "Get buddy works (when no in-group)", {
    state_buddy = state
    group = c(1, 2, 3)
    names(group) = nms
    state_buddy$group = group
    ref = m4ma::getBuddy_r(n, group, a, p_pred, centres,
                           objects, FALSE, state_buddy)
    res = m4ma::getBuddy_rcpp(n, group, a, p_pred, centres,
                           objects, FALSE, state_buddy)
    testthat::expect_equal(res, ref)
  }
)

testthat::test_that(
  "Get buddy returns NULL (when no in-group)", {
    state_buddy = state
    group = c(1, 2, 3)
    names(group) = nms
    state_buddy$group = group
    res = m4ma::getBuddy_rcpp(n, group, a, p_pred, centres,
                              objects, FALSE, state_buddy)
    testthat::expect_equal(res, NULL)
  }
)

testthat::test_that(
  "Get buddy works (pick best)", {
    state_buddy = state
    group = c(2, 2, 2)
    names(group) = nms
    state_buddy$group = group
    p2_same = rbind(c(0, 0), c(0.5, 0), c(0, 0.5))
    rownames(p2_same) = nms
    state_buddy$p = p2_same
    ref = m4ma::getBuddy_r(n, group, a, p_pred, centres,
                           objects, TRUE, state_buddy)
    res = m4ma::getBuddy_rcpp(n, group, a, p_pred, centres,
                           objects, TRUE, state_buddy)
    testthat::expect_equal(res, ref)
  }
)
