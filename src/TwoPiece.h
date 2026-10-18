
#ifndef GAS_TWO_PIECE_H
#define GAS_TWO_PIECE_H
#include <RcppArmadillo.h>
#include <cmath>
struct GAS_TwoPiece {
  double mean, sd, logleft, logright;
  GAS_TwoPiece(double xi, double m1) {
    const double delta = xi - 1.0/xi;
    mean = m1*delta;
    sd = sqrt(1.0 + (1.0-m1*m1)*delta*delta);
    const double a = 2.0*log(xi);
    const double sp = a > 0 ? a + log1p(exp(-a)) : log1p(exp(a));
    logleft = -sp;
    logright = a-sp;
  }
};
inline double GAS_t_absmean(double nu) {
  return exp(log(2.0) + 0.5*log(nu-2.0) - log(nu-1.0)
    - Rf_lbeta(0.5, nu/2.0));
}
#endif
