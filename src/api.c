#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include "AUCt.h"
#include "nrutil.h"

SEXP survivalauc_compute_auc_components(SEXP time,
                                        SEXP status,
                                        SEXP score,
                                        SEXP threshold_in,
                                        SEXP unique_time) {
  // Replace with call to `AUCt()` and massaging of inputs/results

  double *Y, *D, *Z, *uniqY, *threshold;
  Y = REAL(time);
  D = REAL(status);
  Z = REAL(score);
  threshold = REAL(threshold_in);
  uniqY = REAL(unique_time);

  int n, m, k;
  n = Rf_length(time);
  m = Rf_length(unique_time);
  k = Rf_length(threshold_in);

  double **sens, **spec, *auc;
  sens = dmatrix(1, n, 1, k);
  spec = dmatrix(1, n, 1, k);
  auc = dvector(1, n);

  AUCt(&n, Y, D, Z, &k, threshold, &m, uniqY, sens, spec, auc);

  return R_NilValue;
}
