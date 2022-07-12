#' Probability of the Conditional Nested Logit Model
#'
#' @param cell Alternative vector with nest indices.
#' @param V Vector with utility for each alternative.
#' @param muM Vector with nests associations ranging between 0 and 1.
#' @param nests List of vectors with utility indices.
#' @param alpha List of vectors with alpha values.
#' @param mu General nest association ranging between 0 and 1.
#' @param cellNest Integer matrix associating the cells with the nests.
#'
#' @return Probability of alternative \code{cell} given \code{V}, \code{alpha}, and parameters \code{muM} and \code{mu}.
#' @export
pCNL <- function(cell, V, muM = rep(1, length(nests)), nests, alpha, mu, cellNest) {
  
  # Probability of alternatives within nests
  pAinNest <- function(Vlist, nests, alpha, muM) {
    pim <- nests
    for (m in 1:length(nests)) { # alternatives in each nest
      tmp <- alpha[[m]] * exp(muM[m] * Vlist[[m]])  
      bad <- tmp == Inf # Overflow
      if (any(bad)) tmp <- ifelse(bad, 1, 0)
      if (all(tmp == 0)) {
        pim[[m]] <- tmp 
      } else {
        pim[[m]] <- tmp / sum(tmp)
      }
    }
    return(pim)
  }
  
  # Probability of nest m
  pNest <- function(Vlist, nests, alpha, mu, muM) {
    mu_muM <- mu / muM
    tmp <- sapply(1:length(nests), function(m) {
      (sum(alpha[[m]] * exp(muM[m] * Vlist[[m]])))^mu_muM[m]
    })
    bad <- tmp == Inf # Overflow
    if (any(bad)) {
      tmp <- ifelse(bad, 1, 0)
    }
    if (all(tmp == 0)) {
      tmp 
    } else {
      tmp / sum(tmp)
    }
  }
  
  # Add 1 to cells if there is a cell for stopping, 0 can't be used for index
  if (any(unlist(nests) == 0)) {
    nests <- lapply(nests, function(x) {
      x + 1
    })
    cell <- cell + 1
  }
  
  # Nest list with utility of cells in nest
  Vlist <- lapply(nests, function(x) {
    V[x]
  })
  
  # Set largest V to zero to avoid numerical issues
  maxV <- max(unlist(lapply(Vlist, max)))
  Vlist <- lapply(Vlist, function(x) {
    x - maxV
  })
  
  # Probability of nest
  pN <- pNest(Vlist, nests, alpha, mu, muM)
  
  # Probability of alternative in nest
  pAN <- pAinNest(Vlist, nests, alpha, muM)
  
  # pAN * pN [NonCentral/Central] + pAN * pN [acc/const/dec]
  pAN[[cellNest[[cell, 1]]]][cellNest[cell, 3]] *    # [[ind nest]][ind cell]
    pN[cellNest[cell, 1]] +                          # [ind nest]
    pAN[[cellNest[[cell, 2]]]][cellNest[cell, 4]] *  # [[ind nest]][ind cell]
    pN[cellNest[cell, 2]]                            # [ind nest]
}
