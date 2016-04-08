#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat invSympdC(arma::mat m) {
	return arma::inv_sympd(m);
}

// [[Rcpp::export]]
arma::mat invC(arma::mat m) {
	return arma::inv(m);
}

// [[Rcpp::export]]
arma::mat symC(arma::mat m) {
	return (m + m.t()) / 2;
}
