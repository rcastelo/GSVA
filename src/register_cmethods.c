#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

/* prototypes of functions to be registered */

SEXP
matrix_density_R(double* X, double* Y, int* n_density_samples,
                 int* n_test_samples, int* n_genes, int* rnaseq);

SEXP
ks_matrix_R(SEXP XR, SEXP sidxsR, SEXP n_genesR, SEXP geneset_idxsR,
            SEXP n_genesetR, SEXP tauR, SEXP n_samplesR, SEXP mx_diffR, SEXP abs_rnkR);

/* registration of C-entry points */

static R_CallMethodDef callMethods[] = {
  {"ks_matrix_R", (DL_FUNC) &ks_matrix_R, 9},
  {"matrix_density_R", (DL_FUNC) &matrix_density_R, 6},
  {NULL, NULL, 0}
};

void
R_init_GSVA(DllInfo *info) {

  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);

}
