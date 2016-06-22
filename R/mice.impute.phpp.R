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
#' Garland Jr, Theodore, and Anthony R. Ives. Using the past to predict
#' the present: confidence intervals for regression equations in phylogenetic
#' comparative methods. \emph{The American Naturalist}, 155.3 (2000): 346-364.
#' \url{http://www.jstor.org/stable/10.1086/303327}
#' @export
mice.impute.phpp <- function(y, ry, x, psi, ...) {
  x <- cbind(1, as.matrix(x))

  mat2 <- psi/max(psi)
  k <- 2
  mat3 <- k - mat2
  mat4 <- apply(mat3, 1, function (x) x/sum(x))

  vals4 <- mat4[, 4]
  vals44 <- vals4[-4]
  vals44/sum(vals44)
  donortip <- sample(c(1:3, 5:8), 1, prob=vals44/sum(vals44))

  ymiss <- y[!ry]
  C <- psi[!ry, !ry]
  nmiss <- sum(!ry)

}


