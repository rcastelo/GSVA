#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* prototypes of functions to be registered */

void
matrix_density_R(double* X, double* Y, double* R, int* n_density_samples,int* n_test_samples, int* n_genes, int* rnaseq);

void
ks_matrix_R(double* X, double* R, int* sidxs, int* n_genes, int* geneset_idxs, int* n_geneset, double* tau,  int* n_samples, int* mx_diff);

/* registration of C-entry points */

static R_NativePrimitiveArgType
matrix_density_R_t[7] = {REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP};

static R_NativePrimitiveArgType
ks_matrix_R_t[9] = {REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP};

static const R_CMethodDef
cMethods[] = {
  {"matrix_density_R", (DL_FUNC) &matrix_density_R, 7, matrix_density_R_t},
  {"ks_matrix_R", (DL_FUNC) &ks_matrix_R, 9, ks_matrix_R_t},
  {NULL, NULL, 0}
};

void
R_init_GSVA(DllInfo *info) {
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}
