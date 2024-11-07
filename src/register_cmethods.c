#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

/* prototypes of functions to be registered */

SEXP
matrix_density_R(SEXP density_dataR, SEXP test_dataR, SEXP n_density_samplesR,
                 SEXP n_test_samplesR, SEXP n_genesR, SEXP GausskR,
                 SEXP any_naR, SEXP na_useR, SEXP verboseR);

SEXP
ks_matrix_R(SEXP XR, SEXP sidxsR, SEXP n_genesR, SEXP geneset_idxsR,
            SEXP n_genesetR, SEXP tauR, SEXP n_samplesR, SEXP mx_diffR, SEXP abs_rnkR);

SEXP
kcdfvals_sparse_to_sparse_R(SEXP XCspR, SEXP XRspR, SEXP GausskR, SEXP verboseR);

SEXP
kcdfvals_sparse_to_dense_R(SEXP XCspR, SEXP XRspR, SEXP GausskR, SEXP verboseR);

SEXP
ecdfvals_sparse_to_sparse_R(SEXP XCspR, SEXP XRspR, SEXP verboseR);

SEXP
ecdfvals_sparse_to_dense_R(SEXP XCspR, SEXP XRspR, SEXP verboseR);

SEXP
ecdfvals_dense_to_dense_R(SEXP XR, SEXP verboseR);

SEXP
ecdfvals_dense_to_dense_nas_R(SEXP XR, SEXP verboseR);

SEXP
order_rankstat_R(SEXP xR);

SEXP
gsva_score_genesets_R(SEXP genesetsidxR, SEXP decordstatR, SEXP symrnkstatR,
                      SEXP maxdiffR, SEXP absrnkR, SEXP tauR, SEXP anynaR,
                      SEXP nauseR, SEXP minsizeR);

SEXP
order_rankstat_sparse_to_dense_R(SEXP XCspR, SEXP jR);

SEXP
order_rankstat_sparse_to_sparse_R(SEXP XCspR, SEXP jR);

/* registration of C-entry points */

static R_CallMethodDef callMethods[] = {
  {"ks_matrix_R", (DL_FUNC) &ks_matrix_R, 9},
  {"matrix_density_R", (DL_FUNC) &matrix_density_R, 9},
  {"kcdfvals_sparse_to_sparse_R", (DL_FUNC) &kcdfvals_sparse_to_sparse_R, 4},
  {"kcdfvals_sparse_to_dense_R", (DL_FUNC) &kcdfvals_sparse_to_dense_R, 4},
  {"ecdfvals_sparse_to_sparse_R", (DL_FUNC) &ecdfvals_sparse_to_sparse_R, 3},
  {"ecdfvals_sparse_to_dense_R", (DL_FUNC) &ecdfvals_sparse_to_dense_R, 3},
  {"ecdfvals_dense_to_dense_R", (DL_FUNC) &ecdfvals_dense_to_dense_R, 2},
  {"ecdfvals_dense_to_dense_nas_R", (DL_FUNC) &ecdfvals_dense_to_dense_nas_R, 2},
  {"gsva_scores_genesets_R", (DL_FUNC) &gsva_score_genesets_R, 9},
  {"order_rankstat_sparse_to_dense_R", (DL_FUNC) &order_rankstat_sparse_to_dense_R, 2},
  {"order_rankstat_sparse_to_sparse_R", (DL_FUNC) &order_rankstat_sparse_to_sparse_R, 2},
  {NULL, NULL, 0}
};

/* global variables */
SEXP Matrix_DimNamesSym,
     Matrix_DimSym,
     Matrix_xSym,
     Matrix_iSym,
     Matrix_jSym,
     Matrix_pSym;

void
R_init_GSVA(DllInfo *info) {

  R_registerRoutines(info, NULL, callMethods, NULL, NULL);

  /* from the Matrix package init.c */
  Matrix_DimNamesSym = install("Dimnames");
  Matrix_DimSym = install("Dim");
  Matrix_xSym = install("x");
  Matrix_iSym = install("i");
  Matrix_jSym = install("j");
  Matrix_pSym = install("p");

  R_useDynamicSymbols(info, TRUE);

}
