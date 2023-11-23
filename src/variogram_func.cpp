// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
using namespace Rcpp;

//' Compute empirical variogram
//'
//' @param y vector or matrix of values from which to calculate the variogram
//' @param w_mat matrix of weights corresponding to the pairs of dimensions
//' @param h the spatial lag at which to calculate the empirical variogram
//' @param p exponent used when calculating the variogram
//'
//' @details
//' \code{vario()} takes a vector \code{y} of length $d$ and a matrix \code{w_mat}
//' of dimension $(d, d)$. The elements of \code{w_mat} are weights, and must
//' therefore be non-negative.
//'
//' This calculates the double sum from $i = 1, \dots, d$ and $j = 1, \dots, d$
//' of the weighted difference between the components of \code{y}.
//'
//' \code{vario_mat()} takes a matrix \code{y} of dimension $(p, q)$ and a vector
//' \code{h} of length 2. \code{h} represents the spatial lag at which the spatial
//' variogram will be calculated, so the first element of \code{h} must be an integer
//' between 0 and $p$, and the second element must be an integer between 0 and $q$.
//'
//' In both cases, \code{p} is a single positive number that represents the exponent
//' used when calculating the variogram. This is typically set equal to 2.
//'
//' @return
//' The variogram
//'
//' @name variogram_func
//' @importFrom Rcpp evalCpp

//' @rdname variogram_func
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


//' @rdname variogram_func
// [[Rcpp::export]]
double vario_mat(arma::mat y, arma::ivec h, double p){

  int count = 0;
  double out = 0;
  double nr = y.n_rows;
  double nc = y.n_cols;
  for (int i = 0; i < nr; i++) {
    int ih = i + h[0];
    for (int j = 0; j < nc; j++) {
      int jh = j + h[1];
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
