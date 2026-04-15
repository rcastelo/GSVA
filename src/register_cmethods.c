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
kcdfvals_sparse_to_sparse_R(SEXP XCspR, SEXP XRspR, SEXP GausskR, SEXP verboseR);

SEXP
kcdfvals_sparse_to_dense_R(SEXP XCspR, SEXP XRspR, SEXP GausskR, SEXP verboseR);

SEXP
kcdfvals_svt_to_dense_R(SEXP XsvtR, SEXP GausskR, SEXP verboseR);

SEXP
kcdfvals_svt_to_svt_R(SEXP XsvtR, SEXP GausskR, SEXP verboseR);

SEXP
ecdfvals_sparse_to_sparse_R(SEXP XCspR, SEXP XRspR, SEXP verboseR);

SEXP
ecdfvals_sparse_to_dense_R(SEXP XCspR, SEXP XRspR, SEXP verboseR);

SEXP
ecdfvals_dense_to_dense_R(SEXP XR, SEXP verboseR);

SEXP
ecdfvals_dense_to_dense_nas_R(SEXP XR, SEXP verboseR);

SEXP
ecdfvals_svt_to_dense_R(SEXP XsvtR, SEXP verboseR);

SEXP
ecdfvals_svt_to_svt_R(SEXP XsvtR, SEXP verboseR);

SEXP
gsva_score_genesets_R(SEXP ranksR, SEXP genesetsidxR, SEXP sparseR,
                      SEXP maxdiffR, SEXP absrnkR, SEXP tauR, SEXP anynaR,
                      SEXP nauseR, SEXP minsizeR, SEXP verboseR);

SEXP
ranks2stats_R(SEXP ranksR, SEXP jR, SEXP sparseR);

SEXP
fetch_row_nzvals_R(SEXP svtR, SEXP iR, SEXP whimin1R);

SEXP
row_rngs_nzrngs_RsparseMatrix_R(SEXP XRspR, SEXP verboseR);

/*
SEXP
row_rngs_nzrngs_SVT_SparseMatrix_R(SEXP XsvtR, SEXP verboseR);

SEXP
col_rngs_nzrngs_SVT_SparseMatrix_R(SEXP XsvtR, SEXP verboseR);
*/

SEXP
rowbycols_rngs_nzrngs_SVT_SparseMatrix_R(SEXP XsvtR, SEXP verboseR);

/* registration of C-entry points */

static R_CallMethodDef callMethods[] = {
  {"matrix_density_R", (DL_FUNC) &matrix_density_R, 9},
  {"kcdfvals_sparse_to_sparse_R", (DL_FUNC) &kcdfvals_sparse_to_sparse_R, 4},
  {"kcdfvals_sparse_to_dense_R", (DL_FUNC) &kcdfvals_sparse_to_dense_R, 4},
  {"kcdfvals_svt_to_dense_R", (DL_FUNC) &kcdfvals_svt_to_dense_R, 3},
  {"kcdfvals_svt_to_svt_R", (DL_FUNC) &kcdfvals_svt_to_svt_R, 3},
  {"ecdfvals_sparse_to_sparse_R", (DL_FUNC) &ecdfvals_sparse_to_sparse_R, 3},
  {"ecdfvals_sparse_to_dense_R", (DL_FUNC) &ecdfvals_sparse_to_dense_R, 3},
  {"ecdfvals_svt_to_dense_R", (DL_FUNC) &ecdfvals_svt_to_dense_R, 2},
  {"ecdfvals_svt_to_svt_R", (DL_FUNC) &ecdfvals_svt_to_svt_R, 2},
  {"ecdfvals_dense_to_dense_R", (DL_FUNC) &ecdfvals_dense_to_dense_R, 2},
  {"ecdfvals_dense_to_dense_nas_R", (DL_FUNC) &ecdfvals_dense_to_dense_nas_R, 2},
  {"gsva_score_genesets_R", (DL_FUNC) &gsva_score_genesets_R, 10},
  {"ranks2stats_R", (DL_FUNC) &ranks2stats_R, 4},
  {"fetch_row_nzvals_R", (DL_FUNC) &fetch_row_nzvals_R, 3},
  {"row_rngs_nzrngs_RsparseMatrix_R", (DL_FUNC) &row_rngs_nzrngs_RsparseMatrix_R, 2},
  /* {"row_rngs_nzrngs_SVT_SparseMatrix_R", (DL_FUNC) &row_rngs_nzrngs_SVT_SparseMatrix_R, 2}, */
  /* {"col_rngs_nzrngs_SVT_SparseMatrix_R", (DL_FUNC) &col_rngs_nzrngs_SVT_SparseMatrix_R, 2}, */
  {"rowbycols_rngs_nzrngs_SVT_SparseMatrix_R", (DL_FUNC) &rowbycols_rngs_nzrngs_SVT_SparseMatrix_R, 2},
  {NULL, NULL, 0}
};

/* global variables */
SEXP Matrix_DimNamesSym,
     Matrix_DimSym,
     Matrix_xSym,
     Matrix_iSym,
     Matrix_jSym,
     Matrix_pSym,
     SVT_SparseArray_typeSym,
     SVT_SparseArray_dimNamesSym,
     SVT_SparseArray_dimSym,
     SVT_SparseArray_svtSym,
     GSVA_attrNAsSym;

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

  SVT_SparseArray_typeSym = install("type");
  SVT_SparseArray_dimNamesSym = install("dimnames");
  SVT_SparseArray_dimSym = install("dim");
  SVT_SparseArray_svtSym = install("SVT");

  GSVA_attrNAsSym = install("attrNAs");

  R_useDynamicSymbols(info, TRUE);
}
