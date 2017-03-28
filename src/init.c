// RegisteringDynamic Symbols
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls
*/

extern SEXP GAS_mWCRPS_backtest(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GAS_EvaluateLogScore_Univ(SEXP, SEXP, SEXP, SEXP);
extern SEXP GAS_EvaluateLogScore_Multi(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GAS_ddist_univ(SEXP, SEXP, SEXP, SEXP);
extern SEXP GAS_ddist_multi(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GAS_rdist_univ(SEXP, SEXP);
extern SEXP GAS_rdist_multi(SEXP, SEXP, SEXP);
extern SEXP GAS_pdist_univ(SEXP, SEXP, SEXP);
extern SEXP GAS_qdist_univ(SEXP, SEXP, SEXP);
extern SEXP GAS_mdist_univ(SEXP, SEXP);
extern SEXP GAS_mdist_multi_mean(SEXP, SEXP, SEXP);
extern SEXP GAS_mdist_multi_cov(SEXP, SEXP, SEXP);
extern SEXP GAS_GASFilter_univ(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GAS_GASFilter_multi(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GAS_uGASMultiForcast(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GAS_mGASMultiForcast(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GAS_IM_univ(SEXP, SEXP);
extern SEXP GAS_Map_Vec(SEXP, SEXP, SEXP);
extern SEXP GAS_unmapVec_C(SEXP, SEXP, SEXP);
extern SEXP GAS_MapParameters_univ(SEXP, SEXP, SEXP);
extern SEXP GAS_UnmapParameters_univ(SEXP, SEXP, SEXP);
extern SEXP GAS_MapR_C(SEXP, SEXP);
extern SEXP GAS_UnMapR_C(SEXP, SEXP);
extern SEXP GAS_MapParameters_multi(SEXP, SEXP, SEXP, SEXP);
extern SEXP GAS_UnmapParameters_multi(SEXP, SEXP, SEXP, SEXP);
extern SEXP GAS_EvalMoments_univ(SEXP, SEXP);
extern SEXP GAS_EvalMoments_multi(SEXP, SEXP, SEXP);
extern SEXP GAS_rmvnorm_mat(SEXP, SEXP, SEXP);
extern SEXP GAS_StaticLLK_Univ(SEXP, SEXP, SEXP, SEXP);
extern SEXP GAS_StaticLLK_Multi(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GAS_EvaluatePit_Univ(SEXP, SEXP, SEXP, SEXP);
extern SEXP GAS_Quantiles(SEXP, SEXP, SEXP);
extern SEXP GAS_Score_univ(SEXP, SEXP, SEXP);
extern SEXP GAS_Score_multi(SEXP, SEXP, SEXP, SEXP);
extern SEXP GAS_SimulateGAS_univ(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GAS_SimulateGAS_multi(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GAS_NumberParameters(SEXP, SEXP);
extern SEXP GAS_build_vR(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"GAS_mWCRPS_backtest",         (DL_FUNC) &GAS_mWCRPS_backtest,             8},
  {"GAS_EvaluateLogScore_Univ",   (DL_FUNC) &GAS_EvaluateLogScore_Univ,       4},
  {"GAS_EvaluateLogScore_Multi",  (DL_FUNC) &GAS_EvaluateLogScore_Multi,      5},
  {"GAS_ddist_univ",              (DL_FUNC) &GAS_ddist_univ,                  4},
  {"GAS_ddist_multi",             (DL_FUNC) &GAS_ddist_multi,                 5},
  {"GAS_rdist_univ",              (DL_FUNC) &GAS_rdist_univ,                  2},
  {"GAS_rdist_multi",             (DL_FUNC) &GAS_rdist_multi,                 3},
  {"GAS_pdist_univ",              (DL_FUNC) &GAS_pdist_univ,                  3},
  {"GAS_qdist_univ",              (DL_FUNC) &GAS_qdist_univ,                  3},
  {"GAS_mdist_univ",              (DL_FUNC) &GAS_mdist_univ,                  2},
  {"GAS_mdist_multi_mean",        (DL_FUNC) &GAS_mdist_multi_mean,            3},
  {"GAS_mdist_multi_cov",         (DL_FUNC) &GAS_mdist_multi_cov,             3},
  {"GAS_GASFilter_univ",          (DL_FUNC) &GAS_GASFilter_univ,              8},
  {"GAS_GASFilter_multi",         (DL_FUNC) &GAS_GASFilter_multi,             9},
  {"GAS_uGASMultiForcast",        (DL_FUNC) &GAS_uGASMultiForcast,           10},
  {"GAS_mGASMultiForcast",        (DL_FUNC) &GAS_mGASMultiForcast,           11},
  {"GAS_IM_univ",                 (DL_FUNC) &GAS_IM_univ,                     2},
  {"GAS_Map_Vec",                 (DL_FUNC) &GAS_Map_Vec,                     3},
  {"GAS_unmapVec_C",              (DL_FUNC) &GAS_unmapVec_C,                  3},
  {"GAS_MapParameters_univ",      (DL_FUNC) &GAS_MapParameters_univ,          3},
  {"GAS_UnmapParameters_univ",    (DL_FUNC) &GAS_UnmapParameters_univ,        3},
  {"GAS_MapR_C",                  (DL_FUNC) &GAS_MapR_C,                      2},
  {"GAS_UnMapR_C",                (DL_FUNC) &GAS_UnMapR_C,                    2},
  {"GAS_MapParameters_multi",     (DL_FUNC) &GAS_MapParameters_multi,         4},
  {"GAS_UnmapParameters_multi",   (DL_FUNC) &GAS_UnmapParameters_multi,       4},
  {"GAS_EvalMoments_univ",        (DL_FUNC) &GAS_EvalMoments_univ,            2},
  {"GAS_EvalMoments_multi",       (DL_FUNC) &GAS_EvalMoments_multi,           3},
  {"GAS_rmvnorm_mat",             (DL_FUNC) &GAS_rmvnorm_mat,                 3},
  {"GAS_StaticLLK_Univ",          (DL_FUNC) &GAS_StaticLLK_Univ,              4},
  {"GAS_StaticLLK_Multi",         (DL_FUNC) &GAS_StaticLLK_Multi,             5},
  {"GAS_EvaluatePit_Univ",        (DL_FUNC) &GAS_EvaluatePit_Univ,            4},
  {"GAS_Quantiles",               (DL_FUNC) &GAS_Quantiles,                   3},
  {"GAS_Score_univ",              (DL_FUNC) &GAS_Score_univ,                  3},
  {"GAS_Score_multi",             (DL_FUNC) &GAS_Score_multi,                 4},
  {"GAS_SimulateGAS_univ",        (DL_FUNC) &GAS_SimulateGAS_univ,            6},
  {"GAS_SimulateGAS_multi",       (DL_FUNC) &GAS_SimulateGAS_multi,           7},
  {"GAS_NumberParameters",        (DL_FUNC) &GAS_NumberParameters,            2},
  {"GAS_build_vR",                (DL_FUNC) &GAS_build_vR,                    2},
  {NULL, NULL, 0}
};


void R_init_GAS(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

