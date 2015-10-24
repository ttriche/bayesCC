#' Align cluster label vectors or cluster label matrices as needed
#'
#' @param Z1      either a vector or a matrix
#' @param Z2      either a vector or a matrix
#' @param type    (optional) string: are Z1 & Z2 "vec"tors or "mat"rices? (auto)
#' 
#' @return either a vector or a matrix, depending upon the value of "type" 
#'
#' @export 
#' 
alignClusters <- function(Z1, Z2, type=NULL) {

  if (is.null(type)) {
    type <- ifelse(is.matrix(Z1) & is.matrix(Z2), "mat", "vec")
  } else {
    stopifnot(substr(type, 1, 3) %in% c("vec", "mat"))
  }

  switch(substr(type, 1, 3),
         vec=alignVectors(Z1, Z2),
         mat=alignMatrices(Z1, Z2))

}
