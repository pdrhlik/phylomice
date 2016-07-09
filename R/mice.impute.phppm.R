#' @title phppm
#' @param nn Number of nearest neighbors.
#' @export
mice.impute.phppm <- function(y, ry, x, psi, k = 2, nn = 3, ...) {
  nmiss <- sum(!ry)
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
  res <- vector(length = nmiss)
  j <- 1

  for (i in which(!ry)) {
    vals <- probs[, i]
    valsNoMe <- vals[-i]
    scaledVals <- valsNoMe / sum(valsNoMe)

    # indices of nearest tips
    nearIndices <- c()
    # probabilities of nearest tips
    nearProbs <- c()

    # search for the nn closest tips
    for (ii in 1:nn) {
      m <- which.max(scaledVals)
      nearIndices <- c(nearIndices, m)
      nearProbs <- c(nearProbs, scaledVals[m])
      scaledVals[m] <- 0
    }

    nearProbs <- nearProbs / sum(nearProbs)
    donortip <- sample(nearIndices, 1, prob = nearProbs)

    res[j] <- y[donortip]
    j <- j + 1
  }

  return(res)
}
