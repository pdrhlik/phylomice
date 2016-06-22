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
mice.impute.phnorm <- function(y, ry, x, psi, psiinv, ...) {
	x <- cbind(1, as.matrix(x))
	parm <- .phnorm.draw(y, ry, x, psiinv, ...)

	ymiss <- y[!ry]
	C <- psi[!ry, !ry]
	nmiss <- sum(!ry)

	# vectorized version
	Carr <- lapply(1:nmiss, function(i) invSympdC(C[-i, -i]))
	cihmat <- lapply(1:nmiss, function(i) C[i, -i])
	# heights of the trees
	chhvec <- diag(C)
	# matrix of products of cihmat and Carr
	cihCC <- sapply(1:nmiss, function(i) {
		crossprod(cihmat[[i]], Carr[[i]])
	}, USE.NAMES = FALSE)
	# vector of means
	mu <- sapply(1:nmiss, function(i) {
		y <- ymiss[-i] - .Internal(mean(ymiss[-i]))
		return(crossprod(cihCC[,i], y))
	}, USE.NAMES = FALSE)
	# vector of ch
	ch <- sapply(1:nmiss, function(i) {
		chhvec[i] - crossprod(cihCC[,i], cihmat[[i]])
	}, USE.NAMES = FALSE)

	return(x[!ry, ] %*% parm$beta + rnorm(nmiss, mean = mu, sd = parm$sigma * sqrt(ch)))
}

.phnorm.draw <- function(y, ry, x, psiinv, ridge = 1e-05, ...) {
	xobs <- x[ry, ]
	yobs <- y[ry]
	psiobs <- psiinv[ry, ry]
	ncolX <- ncol(x)

	xProdPsi <- crossprod(xobs, psiobs)

	xtx <- crossprod(t(xProdPsi), xobs)
	pen <- ridge * diag(xtx)
	if (length(pen) == 1) {
		pen <- matrix(pen)
	}

	v <- invC(xtx + diag(pen))
	coef <- crossprod(t(crossprod(v, xProdPsi)), yobs)

	residuals <- yobs - crossprod(t(xobs), coef)
	df <- max(sum(ry) - ncolX, 1)  # SvB 31/10/2012
	sigma.star <- sqrt(sum((residuals)^2)/rchisq(1, df))  # SvB 01/02/2011
	beta.star <- coef + (crossprod(chol(sym(v)), rnorm(ncolX))) * sigma.star
	parm <- list(coef, beta.star, sigma.star)  # SvB 10/2/2010
	names(parm) <- c("coef", "beta", "sigma")  # SvB 10/2/2010
	return(parm)
}

# from internal.R
sym <- function(x) {(x + t(x)) / 2}
