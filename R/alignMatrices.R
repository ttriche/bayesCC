#' Align matrices of cluster labels.
#' 
#' @param Z1  A matrix. 
#' @param Z2  A matrix of the same dimensions as Z1
#' 
#' @return    A matrix of the same dimensions as Z2.
#' 
#' @export 
#' 
alignMatrices <- function(Z1, Z2) { 
  for (k in 1:dim(Z1)[2]) {
    for (tempk in 1:dim(Z2)[2]) {
      Max <- sum(Z1 == Z2)
      Z2dummy <- Z2
      Z2dummy[, k] <- Z2[, tempk]
      Z2dummy[, tempk] <- Z2[, k]
      if (sum(Z1 == Z2dummy) > Max) 
        Z2 <- Z2dummy
    }
  }
  return(Z2) 
}
