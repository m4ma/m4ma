#' c_vd 
#' Calculate cell centres for set of cells (index 1..33) for p1 heading at velocity v1 at angle a1 given time step tStep seconds.
#' @param cells Integer vector indicating the set of pedestrian's cells from 1 to 33
#' @param p1 Numeric vector of pedestrian's xy coordinates 
#' @param v1 Scalar value of current velocity
#' @param a1 current angle
#' @param vels Numeric matrix of velocities
#' @param angles Numeric matrix (33x3) of possible directions expressed in grades from 0 to 360
#' @return Numeric matrix (33x2) of xy coordinates for each cell 
#' @export
c_vd <- function(cells, p1, v1, a1,
                 vels = matrix(rep(c(1.5, 1, .5), each = 11), ncol = 3),
                 angles = matrix(rep(c(72.5, 50, 32.5, 20, 10, 0, 350, 340, 
                                       327.5, 310, 287.5), times = 3), 
                                 ncol = 3)) {
  t(p1 + t(scaleVel(v1) * vels[cells] * aTOd((angles[cells] + a1) %% 360)))
}