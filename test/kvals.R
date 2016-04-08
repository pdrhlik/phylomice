library(ape)
library(phytools)
library(mice)
library(parmice)
library(phylomice)

go <- function(dims, perc, nsteps) {
  tree <- generateTree(dims[1])
  psi <- vcv(tree)
  psiinv <- chol2inv(chol(psi))
  data <- generateData(dims, tree)
  kBefore <- getKValues(data, tree)

  res <- sapply(1:nsteps, function(i) {
    missdata <- generateMissing(data, perc)
    step(missdata, tree, psi, psiinv)
  })
  return(list(
    "kBefore" = kBefore,
    "kAfterMean" = rowMeans(res),
    "kBeforeAfter" = kBefore / rowMeans(res),
    "kAfterAll" = res
  ))
}

step <- function(data, tree, psi, psiinv) {
  imp <- parmice(data, method = "phnorm", ncores = 5, maxit = 100, psi = psi, psiinv = psiinv)
  avgResData <- avgMiceRes(imp)
  kVals <- getKValues(avgResData, tree)
  return(kVals)
}

avgMiceRes <- function(imp) {
  all <- matrix(0, nrow = nrow(imp$data), ncol = ncol(imp$data))
  for (i in 1:imp$m) {
    all <- all + complete(imp, i)
  }
  return(all / imp$m)
}

generateTree <- function(n) {
  return(rtree(n))
}

generateData <- function(dims, tree) {
  sapply(1:dims[2], function(i) {
    rTraitCont(tree)
  })
}

getKValues <- function(data, tree) {
  sapply(1:ncol(data), function(i) {
    phylosig(tree, data[, i], method = "K")
  })
}

generateMissing <- function(data, perc) {
  nmiss <- round(nrow(data) * perc)
  sapply(1:ncol(data), function(i) {
    col <- data[, i]
    col[sample(nrow(data), nmiss)] <- NA
    col
  })
}

