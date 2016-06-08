library(ape)
library(phytools)
library(mice)
library(parmice)
library(phylomice)

runSimulations <- function() {
  setup <- list(
    treeSize = c(32, 64, 128, 512, 1024),
    traitSize = c(2, 5, 10, 15),
    treeShape = c("balanced", "left", "pureBirth"),
    lambdaTransform = c(0, 0.5, 0.8, 1),
    missDataMech = c("mcar", "mar"),
    missDataLevel = c(1, 2, 3),
    miceMethod = c("norm", "phnorm", "pmm", "phpmm")
  )

  setup <- list(
    treeSize = c(32),
    traitSize = c(2),
    treeShape = c("left"),
    lambdaTransform = c(0.5),
    missDataMech = c("mar"),
    missDataLevel = c(1),
    miceMethod = c("phnorm")
  )

  completeResult <- list()

  for (treeSize in setup$treeSize) {
    for (traitSize in setup$traitSize) {
      errorFileName <- paste("error-", traitSize, ".csv", sep = "")
      warningFileName <- paste("warning-", traitSize, ".csv", sep = "")
      resultFileName <- paste("result-", traitSize, ".csv", sep = "")

      writeErrorHeader <- TRUE
      writeWarningHeader <- TRUE
      writeResultHeader <- TRUE
      for (treeShape in setup$treeShape) {
        for (lambdaTransform in setup$lambdaTransform) {
          for (missDataMech in setup$missDataMech) {
            for (missDataLevel in setup$missDataLevel) {
              for (miceMethod in setup$miceMethod) {

                result <- list(
                  "setup" = list(
                    treeSize = treeSize,
                    traitSize = traitSize,
                    treeShape = treeShape,
                    lambdaTransform = lambdaTransform,
                    missDataMech = missDataMech,
                    missDataLevel = missDataLevel,
                    miceMethod = miceMethod
                  )
                )

                print(unlist(result$setup))

                tryCatch({
                  startTime <- Sys.time()
                  tree <- generateTree(treeSize, treeShape)
                  data <- generateData(c(treeSize, traitSize), tree, lambdaTransform)

                  psi <- NULL
                  psiinv <- NULL

                  if (miceMethod == "phnorm" || miceMethod == "phpmm") {
                    psi <- vcv(tree, corr = TRUE)
                    psiinv <- chol2inv(chol(psi))
                  }

                  anBefore <- analyseData(data, tree)
                  result$before = anBefore

                  missdata <- generateMissingByMech(data, missDataMech, missDataLevel)
                  result$numberOfNA <- numberOfNA(missdata)
                  result$numberOfNARows <- numberOfNARows(missdata)

                  startImpTime <- Sys.time()
                  imp <- mice(missdata,
                              method = miceMethod,
                              ncores = 5,
                              maxit = 100,
                              psi = psi,
                              psiinv = psiinv,
                              printFlag = FALSE,
                              version = "2.21")
                  impDuration <- Sys.time() - startImpTime
                  units(impDuration) <- "secs"
                  result$impTime <- impDuration

                  result$imp1 = analyseData(complete(imp, 1), tree)
                  result$imp2 = analyseData(complete(imp, 2), tree)
                  result$imp3 = analyseData(complete(imp, 3), tree)
                  result$imp4 = analyseData(complete(imp, 4), tree)
                  result$imp5 = analyseData(complete(imp, 5), tree)
                }, warning = function(w) {
                  print("WARNING")
                  print(w)
                  result$warning <<- w
                  duration <- Sys.time() - startTime
                  print("===================================")
                  units(duration) <- "secs"
                  result$time <- duration

                  rowRes <- unlist(result, recursive = TRUE)

                  if (writeWarningHeader) {
                    write(paste(names(rowRes), collapse = ";;"), file = warningFileName, append = TRUE, sep = ",")
                    writeWarningHeader <<- FALSE
                  }
                  write(paste(unname(rowRes), collapse = ";;"), file = warningFileName, append = TRUE, sep = ",")
                }, error = function(e) {
                  print("ERROR")
                  print(e)
                  result$error <<- e
                  duration <- Sys.time() - startTime
                  print("===================================")
                  units(duration) <- "secs"
                  result$time <- duration

                  rowRes <- unlist(result, recursive = TRUE)

                  if (writeErrorHeader) {
                    write(paste(names(rowRes), collapse = ";;"), file = errorFileName, append = TRUE, sep = ",")
                    writeErrorHeader <<- FALSE
                  }
                  write(paste(unname(rowRes), collapse = ";;"), file = errorFileName, append = TRUE, sep = ",")
                }, finally = {
                  if (is.null(result$error) && is.null(result$warning)) {
                    duration <- Sys.time() - startTime
                    print("OK")
                    units(duration) <- "secs"
                    result$time <- duration

                    rowRes <- unlist(result, recursive = TRUE)

                    if (writeResultHeader) {
                      write(paste(names(rowRes), collapse = ";;"), file = resultFileName, append = TRUE, sep = ",")
                      writeResultHeader <- FALSE
                    }
                    write(paste(unname(rowRes), collapse = ";;"), file = resultFileName, append = TRUE, sep = ",")

                    print("===================================")
                  }
                })
              }
            }
          }
        }
      }
    }
  }

}

numberOfNA <- function(data) {
  apply(data, 2, FUN = function(x) {
    sum(is.na(x))
  })
}

numberOfNARows <- function(data) {
  res <- apply(data, 1, FUN = function(x) {
    ifelse(sum(is.na(x)) == (length(x) - 1), 1, 0)
  })
  return(sum(res))
}

analyseData <- function(data, tree) {
  result <- list(
    mean = vector(mode = "numeric", length = ncol(data)),
    median = vector(mode = "numeric", length = ncol(data)),
    sd = vector(mode = "numeric", length = ncol(data)),
    mode = vector(mode = "numeric", length = ncol(data)),
    K = vector(mode = "numeric", length = ncol(data)),
    lambda = vector(mode = "numeric", length = ncol(data))
  )
  kValues <- getKValues(data, tree, "K")
  lambdaValues <- getKValues(data, tree, "lambda")
  for (i in 1:ncol(data)) {
    col <- data[, i]
    dens <- density(col[!is.na(col)])
    result$mean[i] <- mean(col)
    result$median[i] = median(col)
    result$sd[i] = sd(col)
    result$mode[i] = dens$x[which.max(dens$y)]
    result$K[i] = kValues[i]
    result$lambda[i] = lambdaValues[i]
  }
  return(result)
}

generateMissingByMech <- function(data, missDataMech, missDataLevel) {
  if (missDataMech == "mcar") {
    perc <- 0.1
    if (missDataLevel == 2) {
      perc <- 0.5
    }
    if (missDataLevel == 3) {
      perc <- 0.8
    }
    return(cbind(data[, 1], generateMissing(matrix(data[, 2:ncol(data)]), perc)))
  } else if (missDataMech == "mar") {
    k <- 0.3
    if (missDataLevel == 2) {
      k <- 0.8
    }
    if (missDataLevel == 3) {
      k <- 1.3
    }
    x1 <- data[, 1]
    probs <- exp(-k*abs(x1)) # k is a tuning parameter governing the strength of the relationship with x1
    for (i in 2:ncol(data)) {
      data[, i] <- ifelse(probs < runif(nrow(data)), NA, data[, i])
    }
    return(data)
  }
}

go <- function(dims, perc, nsteps, method = "lambda") {
  cat("SETTING UP DATA...\n")
  start <- Sys.time()
  tree <- generateTree(dims[1])
  psi <- vcv(tree, corr = TRUE)
  psiinv <- chol2inv(chol(psi))
  data <- generateData(dims, tree, 0.8)
  dataNames <- names(data)
  kBefore <- getKValues(data, tree, method)
  missdata <- generateMissing(data, perc)
  str(psi)
  str(psiinv)
  str(missdata)
  str(tree)
  res <- sapply(1:nsteps, function(i) {
    cat(Sys.time() - start, ": RUNNING STEP", i, "...\n")
    step(missdata, tree, psi, psiinv, method)
  })
  cat(Sys.time() - start, ": FINISHED", "\n")
  kAfterMean <- rowMeans(res)
  return(list(
    "kBefore" = kBefore,
    "kAfterMean" = kAfterMean,
    "kBeforeAfter" = kBefore / kAfterMean,
    "kResiduals" = sqrt((kBefore - kAfterMean)^2),
    "kAfterAll" = res,
    "method" = method
  ))
}

step <- function(data, tree, psi, psiinv, method) {
  imp <- mice(data, method = "phnorm", ncores = 5, maxit = 100, psi = psi, psiinv = psiinv, printFlag = FALSE, version = "2.21")
  avgKVals <- avgMiceRes(imp, tree, method)
  # kVals <- getKValues(avgResData, tree)
  return(avgKVals)
}

avgMiceRes <- function(imp, tree, method) {
  # all <- matrix(0, nrow = nrow(imp$data), ncol = ncol(imp$data))
  all <- vector(mode = "numeric", length = ncol(imp$data))
  for (i in 1:imp$m) {
    all <- all + getKValues(complete(imp, i), tree, method)
  }
  return(all / imp$m)
}

generateTree <- function(n, treeShape = "pureBirth") {
  print(treeShape)
  if (treeShape == "pureBirth") {
    return(pbtree(n = n))
  } else if (treeShape == "balanced") {
    tree <- stree(n, "balanced")
    return(compute.brlen(tree))
  } else if (treeShape == "left") {
    tree <- stree(n, "left")
    return(compute.brlen(tree))
  }
}

generateData <- function(dims, tree, lambdaTransform) {
  mat <- vcv(tree)
  trans <- lambda.transform(lambdaTransform, mat)
  newTree <- vcv2phylo(trans)
  sapply(1:dims[2], function(i) {
    fastBM(newTree)
  })
}

getKValues <- function(data, tree, method) {
  dm <- dimnames(data)[[1]]
  res <- sapply(1:ncol(data), function(i) {
    a <- data[, i]
    names(a) <- dm
    phylosig(tree, a, method = method)
  })
  if (method == "K") {
    return(res)
  } else if (method == "lambda") {
    return(unlist(res[1, ]))
  }
}

generateMissing <- function(data, perc) {
  nmiss <- round(nrow(data) * perc)
  sapply(1:ncol(data), function(i) {
    col <- data[, i]
    col[sample(nrow(data), nmiss)] <- NA
    col
  })
}
