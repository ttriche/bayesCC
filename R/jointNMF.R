#' attempt a joint nonnegative matrix factorization of two matrices,
#' as described in Wang et al. 2014, doi:10.1093/bioinformatics/btu679
#' 
#' @param mat1    the first matrix
#' @param mat2    the second matrix
#' @param findK   try to find the rank K for each matrix?
#' @param howNA   for rank finding, add NAs column-wise, row-wise, or both?
#' @param viaCV   for rank finding, should five-fold CV be used when imputing?
#' @param fracNA  for rank finding, what fraction of the data should be NA'ed?
#'
#' @import RcppML
#'
#' @export
jointNMF <- function(mat1, mat2, findK=FALSE, howNA=c("both","column","row"),
                     viaCV=FALSE, fracNA=0.2) {
  stop("see RcppML joint NMF example and please submit a PR")
}
