#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

// .Call calls
extern SEXP survivalauc_compute_auc_components(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"survivalauc_compute_auc_components", (DL_FUNC) &survivalauc_compute_auc_components, 4},
  {NULL, NULL, 0}
};

void R_init_survivalauc(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
