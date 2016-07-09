#' @title phpmm
#' @export
mice.impute.phpmm <- function(y, ry, x, psi, psiinv, donors = 5, type = 1,
                            ridge = 1e-05, version = "", ...)
{
    x <- cbind(1, as.matrix(x))
    ynum <- y
    if (is.factor(y)) ynum <- as.integer(y)  ## added 25/04/2012
    parm <- .phnorm.draw(ynum, ry, x, psiinv, ridge = ridge, ...)  ## bug fix 10apr2013
    yhatobs <- x[ry, ] %*% parm$coef
    yhatmis <- x[!ry, ] %*% parm$beta
    if (version == "2.21")
        return(apply(as.array(yhatmis), 1,
                     .pmm.match,
                     yhat = yhatobs,
                     y = y[ry],
                     donors = donors, ...))
    idx <- matcher(yhatobs, yhatmis, k = donors)
    return(y[ry][idx])
}
