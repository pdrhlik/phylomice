#' @param psi Inverse matrix containing phylogenetic information.
#' This matrix should be passed directly to the mice function.
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
		cihC <- crossprod(cihmat[[i]], Carr[[i]])
		y <- ymiss[-i] - .Internal(mean(ymiss[-i]))
		return(diag(crossprod(cihCC[,i], y)))
	}, USE.NAMES = FALSE)
	# vector of ch
	ch <- sapply(1:nmiss, function(i) {
		chhvec[i] - crossprod(cihCC[,i], cihmat[[i]])
	}, USE.NAMES = FALSE)

# 	# nonvectorized version
# 	mu <- vector(mode = "numeric", length = nmiss)
# 	ch <- vector(mode = "numeric", length = nmiss)
# 	for (i in 1:nmiss) {
# 		chh <- C[i,i]
# 		CCinv <- chol2inv(chol(C[-i, -i]))
# 		yy <- ymiss[-i]
#
# 		cih <- C[i, -i]
# 		mu[i] <- cih %*% CCinv %*% (yy - .Internal(mean(yy)))
# 		ch[i] <- chh - t(cih) %*% CCinv %*% cih
# 	}

	return(x[!ry, ] %*% parm$beta + rnorm(nmiss, mean = mu, sd = parm$sigma * ch))
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
