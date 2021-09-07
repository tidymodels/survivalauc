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

  double *Y, *D, *Z, *uniqY, *threshold;
  // `-1` is added here to make the vectors behave as 1-index to match the
  // expected behavior of numerical recipes
  Y = REAL(time)-1;
  D = REAL(status)-1;
  Z = REAL(score)-1;
  threshold = REAL(threshold_in)-1;
  uniqY = REAL(unique_time)-1;

  int n, m, k;
  n = Rf_length(time);
  m = Rf_length(unique_time);
  k = Rf_length(threshold_in);

  double **sens, **spec, *auc;
  sens = dmatrix(1, n, 1, k);
  spec = dmatrix(1, n, 1, k);
  auc = dvector(1, n);

  AUCt(&n, Y, D, Z, &k, threshold, &m, uniqY, sens, spec, auc);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));

  SEXP out_names = PROTECT(Rf_allocVector(STRSXP, 4));
  SET_STRING_ELT(out_names, 0, Rf_mkChar("unique_times"));
  SET_STRING_ELT(out_names, 1, Rf_mkChar("sensitivity"));
  SET_STRING_ELT(out_names, 2, Rf_mkChar("specificity"));
  SET_STRING_ELT(out_names, 3, Rf_mkChar("auc"));
  Rf_setAttrib(out, R_NamesSymbol, out_names);

  SEXP uniqY_out = Rf_allocVector(REALSXP, m);
  SET_VECTOR_ELT(out, 0, uniqY_out);
  double* v_uniqY_out = REAL(uniqY_out);

  for (int i = 1; i <= m; ++i) {
    v_uniqY_out[i-1] = uniqY[i];
  }

  SEXP sens_out = Rf_allocMatrix(REALSXP, m, k);
  SET_VECTOR_ELT(out, 1, sens_out);
  double* v_sens_out = REAL(sens_out);
  R_xlen_t sens_out_loc = 0;

  for (int col = 1; col <= k; ++col) {
    for (int row = 1; row <= m; ++row) {
      // sens[row=1][col=1], sens[row=2][col=1], ...
      v_sens_out[sens_out_loc] = sens[row][col];
      ++sens_out_loc;
    }
  }

  SEXP spec_out = Rf_allocMatrix(REALSXP, m, k);
  SET_VECTOR_ELT(out, 2, spec_out);
  double* v_spec_out = REAL(spec_out);
  R_xlen_t spec_out_loc = 0;

  for (int col = 1; col <= k; ++col) {
    for (int row = 1; row <= m; ++row) {
      // spec[row=1][col=1], spec[row=2][col=1], ...
      v_spec_out[spec_out_loc] = spec[row][col];
      ++spec_out_loc;
    }
  }

  SEXP auc_out = Rf_allocVector(REALSXP, m);
  SET_VECTOR_ELT(out, 3, auc_out);
  double* v_auc_out = REAL(auc_out);

  for (int i = 1; i <= m; ++i) {
    v_auc_out[i-1] = auc[i];
  }

  UNPROTECT(2);

  return out;
}
