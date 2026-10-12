#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>

/* global variables */
extern SEXP Matrix_DimNamesSym,
            Matrix_DimSym,
            Matrix_xSym,
            Matrix_iSym,
            Matrix_jSym,
            Matrix_pSym,
            SVT_SparseArray_typeSym,
            SVT_SparseArray_dimNamesSym,
            SVT_SparseArray_dimSym,
            SVT_SparseArray_svtSym;

/* accumulate the value 'x' of row 'r' into the statistics 'st' of 'nr' rows */
static inline void
acc_row_stats(double x, int r, double* st, int nr, int logsums) {
  if (ISNAN(x)) {        /* missing value */
    st[nr * 3 + r] += 1;
    return;
  }
  if (x == 0)            /* explicitly stored zero */
    return;
  if (ISNAN(st[r]) || x < st[r])
    st[r] = x;
  if (ISNAN(st[nr + r]) || x > st[nr + r])
    st[nr + r] = x;
  st[nr * 2 + r] += 1;
  if (logsums)
    st[nr * 4 + r] += log(x);
}

/* statistics of the rows of a sparse matrix, either an SVT_SparseMatrix or a
 * dgCMatrix object, calculated by walking through its columns: minimum and
 * maximum nonzero value, number of nonzero non-missing values, number of
 * missing values and, optionally, sum of the logarithm of the nonzero
 * non-missing values. these statistics can be accumulated across blocks of
 * columns, see R function .rowStats() */
SEXP
rowbycols_stats_sparse_R(SEXP XR, SEXP svtR, SEXP logsumsR) {
  SEXP     stR;
  Rboolean svt=asLogical(svtR);
  int      logsums=asLogical(logsumsR);
  int*     dim;
  int      nr, nc;
  double*  st;

  if (svt)
    dim = INTEGER(GET_SLOT(XR, SVT_SparseArray_dimSym));
  else
    dim = INTEGER(GET_SLOT(XR, Matrix_DimSym));
  nr = dim[0];
  nc = dim[1];

  PROTECT(stR = allocMatrix(REALSXP, nr, logsums ? 5 : 4));
  st = REAL(stR);
  for (int r=0; r < nr; r++) {
    st[r] = st[nr + r] = NA_REAL;  /* nonzero minimum and maximum */
    st[nr * 2 + r] = st[nr * 3 + r] = 0;
    if (logsums)
      st[nr * 4 + r] = 0;
  }

  if (svt) {
    SEXP Xsvt_SVT = GET_SLOT(XR, SVT_SparseArray_svtSym);
    const char* type = CHAR(STRING_ELT(getAttrib(XR, SVT_SparseArray_typeSym), 0));
    int itypevals = strcmp(type, "double") != 0;

    if (length(Xsvt_SVT) > 0) /* otherwise the input matrix is empty */
      for (int j=0; j < nc; j++) {
        SEXP svtLeaf = VECTOR_ELT(Xsvt_SVT, j);

        if (svtLeaf == R_NilValue)
          continue;

        SEXP valsR = VECTOR_ELT(svtLeaf, 0);
        SEXP offsetsR = VECTOR_ELT(svtLeaf, 1);
        int* offsets = INTEGER(offsetsR);
        int  noffsets = length(offsetsR);

        if (valsR == R_NilValue) /* lacunar leaf, all nonzero values are 1 */
          for (int k=0; k < noffsets; k++)
            acc_row_stats(1.0, offsets[k], st, nr, logsums);
        else if (itypevals) {
          int* vals = INTEGER(valsR);
          for (int k=0; k < noffsets; k++)
            acc_row_stats(vals[k] == NA_INTEGER ? NA_REAL : (double) vals[k],
                          offsets[k], st, nr, logsums);
        } else {
          double* vals = REAL(valsR);
          for (int k=0; k < noffsets; k++)
            acc_row_stats(vals[k], offsets[k], st, nr, logsums);
        }
      }
  } else {
    int*    p = INTEGER(GET_SLOT(XR, Matrix_pSym));
    int*    i = INTEGER(GET_SLOT(XR, Matrix_iSym));
    double* x = REAL(GET_SLOT(XR, Matrix_xSym));

    for (int j=0; j < nc; j++)
      for (int k=p[j]; k < p[j+1]; k++)
        acc_row_stats(x[k], i[k], st, nr, logsums);
  }

  UNPROTECT(1); /* stR */

  return(stR);
}
