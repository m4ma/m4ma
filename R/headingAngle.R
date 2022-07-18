# Absolute angular difference between 11 directions with zero cone having
# angle a1 and a2
headingAngle <- function(a2, a1, angles = c(72.5, 50, 32.5, 20, 10, 0, 350, 
                                            340, 327.5, 310, 287.5)) {
  sapply((angles + a1) %% 360, minAngle, a2 = a2)  
}

