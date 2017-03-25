// RegisteringDynamic Symbols
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* .Call calls

extern SEXP MapR_C(SEXP, int);
extern SEXP UnMapR_C(SEXP, int);
extern SEXP rmvt_mat(int, SEXP, SEXP, double);
extern SEXP Quantiles(SEXP, SEXP, SEXP);
extern SEXP NumberParameters(SEXP, int);
extern SEXP build_vR(SEXP, int);

static const R_CallMethodDef CallEntries[] = {
  {"GAS_MapR_C",           (DL_FUNC) &GAS_MapR_C,           2},
  {"GAS_UnMapR_C",         (DL_FUNC) &GAS_UnMapR_C,         2},
  {"GAS_rmvt_mat",         (DL_FUNC) &GAS_rmvt_mat,         4},
  {"GAS_Quantiles",        (DL_FUNC) &GAS_Quantiles,        3},
  {"GAS_NumberParameters", (DL_FUNC) &GAS_NumberParameters, 2},
  {"GAS_build_vR",         (DL_FUNC) &GAS_build_vR,         2},

  {NULL, NULL, 0}
};
 */

void R_init_GAS(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

