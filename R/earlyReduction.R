#' do dimension reduction (via NMF or SVD) before Bayesian consensus clustering
#'
#' if NMF, the rank can be estimated by 5xCV on NAs, though this can be slow.
#' the underlying rationale is that whatever rank K best recovers artificially 
#' missing data (knocked out column-wise, row-wise, or randomly across both) 
#' is the best estimable rank we are likely to recover.  In order to stabilize
#' the estimate of K, we can run 5x cross-validation and rotate the NAs (set at
#' a default of 20% of the entries to facilitate sampling without replacement). 
#'
#' joint NMF can also be requested (as in Wang et al., Bioinformatics 2015, 
#' doi: 10.1093/bioinformatics/btu679) but in this case the ranks can only be 
#' estimated marginally.  Joint rank estimation (and, by extension, optimal 
#' joint imputation for linked views) is an open research topic as best as we 
#' can tell.  if anyone wants to send a patch we will gladly apply it and a 
#' great many people will probably start using it thereafter. 
#' 
#' if SVD, the rank will be whatever the data supports (i.e. min(nrow, ncol)).  
#' 
#' @param mat     a matrix to decompose (columns are samples, rows are features)
#' @param how     one of "NMF" or "SVD"; SVD is likely to be much faster
#' @param mat2    a 2nd matrix to reduce (optional; for joint factorization)
#' @param joint   if using NMF, should joint factorization be attempted?
#' @param findK   if using marginal NMF, should the optimal rank(s) be sought?
#' @param howNA   for rank finding, add NAs column-wise, row-wise, or both?
#' @param viaCV   for rank finding, should five-fold CV be used when imputing?
#' @param fracNA  for rank finding, what fraction of the data should be NA'ed?
#'
#' @return        a list with W, H, and K for each matrix if using NMF,
#'                or a list with D, U, and V for each matrix if using SVD.
#' 
#' @import        NNLM 
#'
#' @export
earlyReduction <- function(mat, how=c("NMF","SVD"), mat2=NULL, joint=FALSE,
                           findK=FALSE, howNA=c("both","column","row"), 
                           viaCV=FALSE, pctNA=0.2) {
  how <- match.arg(how) 
  stop("not done yet")
}
