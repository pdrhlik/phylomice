# # # Prepare test data
# # rm(list=objects())
# # data <- matrix(rnorm(1200), ncol=4)
# #
# # n <- nrow(data)
# # # Add missing data
# # data[sample(n, 30), 2] <- NA
# # data[sample(n, 50), 3] <- NA
# # data[sample(n, 10), 4] <- NA
# #
# # # load mice and ape libraries and phnorm method
# # library(mice)
# # library(ape)
# # library(Rcpp)
# # library(parallel)
# #
# # Rcpp::sourceCpp('testcpp.cpp')
# #
# # source("mice.impute.phnorm.r")
# #
# # # create phylo matrix and its inverse
# # psi <- vcv(rtree(nrow(data)))
# # psiinv <- solve(psi)
# #
# # # mice settings
# # maxit = 100
# #
# # Rprof()
# #
# # # impute with phylo
library(mice)
library(parmice)
library(phylomice)
system.time({
  imp <- mice(data, method = c("", "norm", "norm", "norm"), maxit = maxit, seed = 42, printFlag = FALSE)
})

system.time({
	impphylo <- mice(data, method = c("", "phnorm", "phnorm", "phnorm"),
		psi = psi, psiinv = psiinv, maxit = maxit, printFlag = FALSE, seed = 42)
})

system.time({
  impphylopar <- parmice(data, ncores = 5, method = c("", "phnorm", "phnorm", "phnorm"),
     psi = psi, psiinv = psiinv, maxit = maxit, printFlag = FALSE, seed = 42)
})
# #
# # Rprof(NULL)
# #
# #
# # # impute without phylo
# system.time({
# 	imp <- mice(data, method = c("", "norm", "norm", "norm"), maxit = maxit, seed = 42, printFlag = FALSE)
# })
# #
# # summaryRprof()
