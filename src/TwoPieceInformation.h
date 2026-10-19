#ifndef GAS_TWO_PIECE_INFORMATION_H
#define GAS_TWO_PIECE_INFORMATION_H

#include "TwoPiece.h"

// Expected Fisher information for the Fernandez--Steel two-piece family.
// First use raw coordinates (join, log base scale, log skew [, nu]), then
// transform to the package's (mean, standard deviation, skew [, nu]).
// The Student base here is an ordinary, not variance-standardized, t_nu.
inline arma::mat GAS_two_piece_information(const arma::vec& theta,
                                            bool student) {
  const unsigned int k = student ? 4u : 3u;
  if (theta.n_elem != k || !theta.is_finite() || theta(1) <= 0.0 ||
      theta(2) <= 0.0 || (student && theta(3) <= 2.0))
    Rcpp::stop("Invalid parameters for two-piece Fisher information.");

  const double sigma = theta(1), xi = theta(2);
  const double nu = student ? theta(3) : 0.0;
  const double m1 = student ? GAS_t_absmean(nu) : sqrt(2.0/M_PI);
  const double delta = xi - 1.0/xi;
  const double delta_xi = 1.0 + 1.0/(xi*xi);
  const GAS_TwoPiece tp(xi, m1);
  const double c = tp.sd, m = tp.mean;
  const double c_xi = (1.0-m1*m1)*delta*delta_xi/c;
  const double m_xi = m1*delta_xi;
  const double b = sigma*(student ? sqrt((nu-2.0)/nu) : 1.0)/c;
  const double sum = xi + 1.0/xi;
  const double d = delta/sum;

  arma::mat raw(k,k,arma::fill::zeros);
  if (student) {
    const double density0 = exp(Rf_lgammafn((nu+1.0)/2.0) -
      Rf_lgammafn(nu/2.0) - 0.5*log(M_PI*nu));
    raw(0,0) = (nu+1.0)/((nu+3.0)*b*b);
    raw(0,2) = raw(2,0) = 8.0*density0*(nu+1.0)/((nu+3.0)*sum*b);
    raw(1,1) = 2.0*nu/(nu+3.0);
    raw(1,3) = raw(3,1) = -2.0/((nu+1.0)*(nu+3.0));
    raw(2,3) = raw(3,2) = d*raw(1,3);
    raw(3,3) = 0.25*(Rf_trigamma(nu/2.0)-Rf_trigamma((nu+1.0)/2.0)) -
      (nu+5.0)/(2.0*nu*(nu+1.0)*(nu+3.0));
  } else {
    raw(0,0) = 1.0/(b*b);
    raw(0,2) = raw(2,0) = 4.0*m1/(sum*b);
    raw(1,1) = 2.0;
  }
  raw(1,2) = raw(2,1) = d*raw(1,1);
  raw(2,2) = raw(1,1) + 1.0 - d*d;

  arma::mat H(k,k,arma::fill::zeros);
  H(0,0) = 1.0;
  H(0,1) = -m/c;
  H(0,2) = -sigma*(m_xi/c - m*c_xi/(c*c));
  H(1,1) = 1.0/sigma;
  H(1,2) = -c_xi/c;
  H(2,2) = 1.0/xi;
  if (student) {
    const double m1_nu = m1*(0.5/(nu-2.0) - 1.0/(nu-1.0) +
      0.5*(Rf_digamma((nu+1.0)/2.0)-Rf_digamma(nu/2.0)));
    const double m_nu = m1_nu*delta;
    const double c_nu = -m1*m1_nu*delta*delta/c;
    H(0,3) = -sigma*(m_nu/c - m*c_nu/(c*c));
    H(1,3) = 1.0/(nu*(nu-2.0)) - c_nu/c;
    H(3,3) = 1.0;
  }
  const arma::mat information = H.t()*raw*H;
  return 0.5*(information+information.t());
}

#endif
