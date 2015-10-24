#' convenience function
#'
#' @param l   number, vector, or matrix to be log-sum-exp'ed
#'
#' @return    a numeric value
#' 
#' @export 
#' 
logSum <- function(l) max(l) + log(sum(exp(l - max(l))))
