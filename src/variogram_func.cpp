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

// [[Rcpp::export]]
double vario_mat(arma::mat y, arma::ivec h, double p){

  int count = 0;
  double out = 0;
  double nr = y.n_rows;
  double nc = y.n_cols;
  for (int i = 0; i < nr; i++) {
    int ih = i + h[1];
    for (int j = 0; j < nc; j++) {
      int jh = j + h[0];
      if ((ih < nr) && (ih >= 0) && (jh < nc) && (jh >= 0)) {
        double vy = pow(abs(y(i, j) - y(ih, jh)), p);
        count += 1;
        out += vy;
      }
    }
  }
  double sc_out = out/(2*count);
  return (sc_out);
}
