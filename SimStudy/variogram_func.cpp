// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
using namespace Rcpp;     

// [[Rcpp::export]]
double vsC(arma::colvec y, double p){
  
  double out = 0;
  double d = y.size();
  for (int i = 1; i < (d+1); i++) {
    for (int j = i; j < (d+1); j++) {
      double vy = pow(abs(y[i-1]-y[j-1]), p);
      out += vy;
    }
  }
  return (out);
}