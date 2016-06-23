#' Imputation by GLM and phylogenetic information
#'
#' Imputes univariate continuous missing data using the generalized least square approach.
#' @aliases mice.impute.phnorm phnorm
#' @param y Incomplete data vector of length \code{n}
#' @param ry Vector of missing data pattern (\code{FALSE}=missing,
#' \code{TRUE}=observed)
#' @param x Matrix (\code{n} x \code{p}) of complete covariates.
#' @param psi Covariance matrix containing phylogenetic information.
#' @param psiinv Inverse of \code{psi}.
#' @param ... Other named arguments.
#' @return A vector of length \code{nmis} with imputations.
#' @note \code{mice.impute.phnorm} is based on the
#' \code{mice.impute.norm} method from \code{mice} package and the GLM
#' approach by Garland and Ives (2000).
#' @author Patrik Drhlik, Simon P. Blomberg, 2016
#' @references Van Buuren, S., Groothuis-Oudshoorn, K. (2011). \code{mice}:
#' Multivariate Imputation by Chained Equations in \code{R}. \emph{Journal of
#' Statistical Software}, \bold{45}(3), 1-67.
#' \url{http://www.jstatsoft.org/v45/i03/}
#'
#' @export
mice.impute.phpp <- function(y, ry, x, psi, k = 2, ...) {
  ymiss <- y[!ry]
  psimiss <- psi[!ry, !ry]
  psiobs <- psi[ry, ry]
  nmiss <- sum(!ry)
  nobs <- sum(ry)
  nall <- length(ry)

  # Scaled distance matrix
  scaled <- psi / max(psi)
  # Similarity matrix
  simil <- k - scaled
  # Transition probabilities matrix
  probs <- apply(simil, 1, function (x) {
    x / sum(x)
  })
  # Resulting vector
  res <- vector(mode = "numeric", length = nmiss)
  j <- 1

  for (i in which(!ry)) {
    vals <- probs[, i]
    valsNoMe <- vals[-i]
    scaledVals <- valsNoMe / sum(valsNoMe)
    # Choose one of the values with a certain probability
    donortip <- sample(seq(1:nall)[-i], 1, prob = scaledVals)
    res[j] <- y[donortip]
    j <- j + 1
  }

  return(res)
}


