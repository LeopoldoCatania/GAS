#ifndef GAS_MOMENT_WARNINGS_H
#define GAS_MOMENT_WARNINGS_H
#include <RcppArmadillo.h>
#include <string>

// Normalize unavailable/numerically non-finite outputs to R's NA, not NaN/Inf.
template <typename T>
inline bool GAS_clean_moments(T& values) {
  bool missing = false;
  for (arma::uword j = 0; j < values.n_elem; ++j) {
    if (!R_FINITE(values[j])) { values[j] = NA_REAL; missing = true; }
  }
  return missing;
}

inline void GAS_warn_moments(const std::string& dist,
                             const std::string& unavailable) {
  if (!unavailable.empty())
    Rcpp::warning("%s: unavailable moments (%s); returning NA for affected entries. "
                  "Check moment existence, parameter validity and numerical range.",
                  dist.c_str(), unavailable.c_str());
}

// Rows are mean, variance, skewness and ordinary kurtosis; columns are dates.
inline void GAS_warn_univ_moments(const arma::mat& values,
                                  const std::string& dist) {
  const char* labels[] = {"mean", "variance", "skewness", "kurtosis"};
  std::string missing;
  for (arma::uword j = 0; j < 4; ++j) {
    if (!values.row(j).is_finite()) {
      if (!missing.empty()) missing += ", ";
      missing += labels[j];
    }
  }
  GAS_warn_moments(dist, missing);
}
#endif
