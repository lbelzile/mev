#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP mev_emplik(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mev_EuclideanWeights(SEXP, SEXP);
extern SEXP mev_gloo2cv(SEXP, SEXP, SEXP, SEXP);
extern SEXP mev_gloo3cv(SEXP, SEXP, SEXP, SEXP);
extern SEXP mev_gloocv(SEXP, SEXP, SEXP, SEXP);
extern SEXP mev_ldirfn(SEXP);
extern SEXP mev_mvrnorm(SEXP, SEXP, SEXP);
extern SEXP mev_mvrnorm_arma(SEXP, SEXP, SEXP);
extern SEXP mev_Pickands_emp(SEXP, SEXP, SEXP);
extern SEXP mev_rbilogspec(SEXP, SEXP);
extern SEXP mev_rbrspec(SEXP, SEXP);
extern SEXP mev_RcppExport_registerCCallable();
extern SEXP mev_rdir(SEXP, SEXP, SEXP);
extern SEXP mev_rdirmixspec(SEXP, SEXP, SEXP, SEXP);
extern SEXP mev_rdirspec(SEXP, SEXP, SEXP, SEXP);
extern SEXP mev_rexstudspec(SEXP, SEXP, SEXP);
extern SEXP mev_rhrspec(SEXP, SEXP);
extern SEXP mev_rlogspec(SEXP, SEXP, SEXP);
extern SEXP mev_rmevA1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mev_rmevA2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mev_rmevasy(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mev_rmevspec_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mev_rneglogspec(SEXP, SEXP, SEXP);
extern SEXP mev_rPbilog(SEXP, SEXP, SEXP);
extern SEXP mev_rPBrownResnick(SEXP, SEXP);
extern SEXP mev_rPdir(SEXP, SEXP, SEXP, SEXP);
extern SEXP mev_rPdirmix(SEXP, SEXP, SEXP, SEXP);
extern SEXP mev_rPexstud(SEXP, SEXP, SEXP);
extern SEXP mev_rPHuslerReiss(SEXP, SEXP);
extern SEXP mev_rPlog(SEXP, SEXP, SEXP);
extern SEXP mev_rPneglog(SEXP, SEXP, SEXP);
extern SEXP mev_rPSmith(SEXP, SEXP, SEXP);
extern SEXP mev_rsmithspec(SEXP, SEXP, SEXP);
extern SEXP mev_Zhang_Stephens(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"mev_emplik",                       (DL_FUNC) &mev_emplik,                       7},
    {"mev_EuclideanWeights",             (DL_FUNC) &mev_EuclideanWeights,             2},
    {"mev_gloo2cv",                      (DL_FUNC) &mev_gloo2cv,                      4},
    {"mev_gloo3cv",                      (DL_FUNC) &mev_gloo3cv,                      4},
    {"mev_gloocv",                       (DL_FUNC) &mev_gloocv,                       4},
    {"mev_ldirfn",                       (DL_FUNC) &mev_ldirfn,                       1},
    {"mev_mvrnorm",                      (DL_FUNC) &mev_mvrnorm,                      3},
    {"mev_mvrnorm_arma",                 (DL_FUNC) &mev_mvrnorm_arma,                 3},
    {"mev_Pickands_emp",                 (DL_FUNC) &mev_Pickands_emp,                 3},
    {"mev_rbilogspec",                   (DL_FUNC) &mev_rbilogspec,                   2},
    {"mev_rbrspec",                      (DL_FUNC) &mev_rbrspec,                      2},
    {"mev_RcppExport_registerCCallable", (DL_FUNC) &mev_RcppExport_registerCCallable, 0},
    {"mev_rdir",                         (DL_FUNC) &mev_rdir,                         3},
    {"mev_rdirmixspec",                  (DL_FUNC) &mev_rdirmixspec,                  4},
    {"mev_rdirspec",                     (DL_FUNC) &mev_rdirspec,                     4},
    {"mev_rexstudspec",                  (DL_FUNC) &mev_rexstudspec,                  3},
    {"mev_rhrspec",                      (DL_FUNC) &mev_rhrspec,                      2},
    {"mev_rlogspec",                     (DL_FUNC) &mev_rlogspec,                     3},
    {"mev_rmevA1",                       (DL_FUNC) &mev_rmevA1,                       6},
    {"mev_rmevA2",                       (DL_FUNC) &mev_rmevA2,                       6},
    {"mev_rmevasy",                      (DL_FUNC) &mev_rmevasy,                      7},
    {"mev_rmevspec_cpp",                 (DL_FUNC) &mev_rmevspec_cpp,                 6},
    {"mev_rneglogspec",                  (DL_FUNC) &mev_rneglogspec,                  3},
    {"mev_rPbilog",                      (DL_FUNC) &mev_rPbilog,                      3},
    {"mev_rPBrownResnick",               (DL_FUNC) &mev_rPBrownResnick,               2},
    {"mev_rPdir",                        (DL_FUNC) &mev_rPdir,                        4},
    {"mev_rPdirmix",                     (DL_FUNC) &mev_rPdirmix,                     4},
    {"mev_rPexstud",                     (DL_FUNC) &mev_rPexstud,                     3},
    {"mev_rPHuslerReiss",                (DL_FUNC) &mev_rPHuslerReiss,                2},
    {"mev_rPlog",                        (DL_FUNC) &mev_rPlog,                        3},
    {"mev_rPneglog",                     (DL_FUNC) &mev_rPneglog,                     3},
    {"mev_rPSmith",                      (DL_FUNC) &mev_rPSmith,                      3},
    {"mev_rsmithspec",                   (DL_FUNC) &mev_rsmithspec,                   3},
    {"mev_Zhang_Stephens",               (DL_FUNC) &mev_Zhang_Stephens,               8},
    {NULL, NULL, 0}
};

void R_init_mev(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

