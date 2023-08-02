// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
double vario(arma::colvec y, arma::mat w_mat, double p){

  double out = 0;
  double d = y.size();
  for (int i = 0; i < d; i++) {
    for (int j = 0; j < d; j++) {
      double vy = pow(abs(y[i] - y[j]), p);
      out += 2*w_mat(i, j)*vy;
    }
  }
  return (out);
}
