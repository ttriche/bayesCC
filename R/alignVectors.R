#' Align vectors of cluster labels.
#' 
#' @param Z1  A vector
#' @param Z2  A vector of the same length as Z1 
#' 
#' @return    A vector of the same length as Z2 
#' 
#' @export 
#' 
alignVectors <- function(Z1, Z2, ...) { 
  for (k in 1:length(unique(Z1))) {
    Max <- sum(Z1 == k & Z2 == k)/(0.01 + sum(Z2 == k) + sum(Z1 == k))
    for (tempk in 1:length(unique(Z2))) {
      if (sum(Z1 == k & 
              Z2 == tempk)/(0.01 + sum(Z2 == tempk) + sum(Z1 == k)) > Max) {
        Max <- sum(Z1 == k & Z2 == tempk) / 
               (0.01 + sum(Z2 == tempk) + sum(Z1 == k))
        dummy <- Z2 == k
        Z2[Z2 == tempk] <- k
        Z2[dummy] <- tempk
      }
    }
  }
  return(Z2) 
} 
