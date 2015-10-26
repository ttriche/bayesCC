#' compute mean adjusted adherence, the recommended figure of merit for bayesCC
#' 
#' @param Results     a list (the output from one completed run of bayesCC(...))
#'
#' @return a vector of mean-adjusted adherences for the run, a.k.a. alphaStar
#'
#' @examples
#' 
#' # see ?bayesCC
#' data(Results)
#' alphaStarDist <- data.frame(lapply(Results, alphaStar))
#' boxplot(alphaStarDist, main="Mean-adjusted adherence by number of clusters")
#'
#' # add a violin plot
#' library(vioplot)
#' cols <- c("tomato","springgreen","darkorange","yellow")
#' for(i in seq_along(alphaStarDist)) {
#'   vioplot(alphaStarDist[[i]], col=cols[i], at=i, add=TRUE)   
#' }
#' 
#' @export
#'
alphaStar <- function(Results) {
  
  requiredFields <- c("AlphaVec", "Cbest") 
  if (!all(requiredFields %in% names(Results))) { # {{{ complain
    message("This function requires output directly from a run of bayesCC().")
    message("You seem to have passed it a list of runs, or something else.")
    message("Please ensure your input contains the following fields:")
    for (field in requiredFields) message("  ", field)
    message("Then try again.")
  } # }}}
  
  # first half of samples are burnin
  cycles <- ncol(Results$AlphaVec)
  burnin <- round(cycles / 2) # discard
  sampled <- seq(burnin + 1, cycles)
 
  # k is implicit
  k <- ncol(Results$Cbest)

  # return the full distribution of AlphaStar
  colMeans((Results$AlphaVec[,sampled] - 1 / k) / (1 - 1 / k))

}
