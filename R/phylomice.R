#' phylomice: Multiple imputation using phylogenetic information
#'
#' The phylomice package provides methods for imputing missing
#' data using phylogenetic information. Currently contained
#' methods:
#'
#' \bold{phnorm}
#' \bold{phpmm} - under development
#'
#' @docType package
#' @name phylomice
#'
#' @useDynLib phylomice
#' @import RcppArmadillo
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
NULL

#' Precomputes psi and psiinv given a tree.
#'
#' @param tree Phylogenetic tree that can be processed
#' by the ape library.
#' @return Returns a list with components psi
#' and psiinv that are needed in phylogenetic
#' imputation methods.
#' @export
precomputePsi <- function(tree) {
  psi <- ape::vcv(tree)
  return(list(
    psi = psi,
    psiinv = chol2inv(chol(psi))
  ))
}
