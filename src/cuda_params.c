#include "cuda_params.h"
#include "cuda/config.h"
#include <R.h>
#include <Rdefines.h>
#include <float.h>
#include <math.h>

SEXP gsva_cuda_max_gset_size_R(void) {
  SEXP res;
  PROTECT(res = allocVector(INTSXP, 1));
  INTEGER(res)[0] = GSVA_MAX_GPU_GSET_SIZE;
  UNPROTECT(1);
  return res;
}

// safety of the GPU walk: the walk is safe when the largest possible weight,
// (G/2)^tau, fits in gsva_float_t
SEXP gsva_cuda_walk_safe_R(SEXP GR, SEXP ksR, SEXP tauR) {
  int G = asInteger(GR);
  double tau = asReal(tauR);
  R_xlen_t nk = XLENGTH(ksR);
  int *ks = INTEGER(ksR);

  const double limit =
    log(sizeof(gsva_float_t) == 4 ? (double) FLT_MAX : DBL_MAX) - 1e-3;

  SEXP res;
  PROTECT(res = allocVector(LGLSXP, nk));
  int *rp = LOGICAL(res);
  for (R_xlen_t i = 0; i < nk; i++) {
    if (!R_finite(tau) || tau < 0.0) {
      rp[i] = FALSE;  // defensive only
    } else if (ks[i] < 1 || G <= 2) {
      rp[i] = TRUE;  // empty set or unit rank statistics cannot saturate
    } else {
      double log_bound = tau * log((double) G / 2.0);
      rp[i] = (log_bound < limit) ? TRUE : FALSE;
    }
  }

  UNPROTECT(1);
  return res;
}
