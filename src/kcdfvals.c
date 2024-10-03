#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <cli/progress.h>

/* from Mutils.h */
static R_INLINE
SEXP ALLOC_SLOT(SEXP obj, SEXP nm, SEXPTYPE type, int length)
{
  SEXP val = allocVector(type, length);

  SET_SLOT(obj, nm, val);
  return val;
}

/* global variables */
extern SEXP Matrix_DimNamesSym,
            Matrix_DimSym,
            Matrix_xSym,
            Matrix_iSym,
            Matrix_jSym,
            Matrix_pSym;

double
sd(double* x, int n);

void
outerselfsubtr(double* x, int n, double *out); 

void
row_d_nologodds(double* x, double* y, double* r, int size_density_n,
                int size_test_n, int rnaseq);


SEXP
kcdfvals_sparse_to_sparse_R(SEXP XCspR, SEXP XRspR, SEXP GausskR, SEXP verboseR);

SEXP
kcdfvals_sparse_to_dense_R(SEXP XCspR, SEXP XRspR, SEXP GausskR, SEXP verboseR);

/*
SEXP
kcdfvals_dense_to_dense_R(SEXP XR);
*/

SEXP
allocMatrix(SEXPTYPE mode, int nrow, int ncol);

/* numerical (double) comparison for qsort() */
extern int
dbl_cmp(const void* a, const void* b);

/* calculate the outer product of a vector with itself
 * using subtraction as the arithmetic operator, only
 * among values in different positions in 'x'. the output
 * is a vector conforming to a column-major lower-triangular
 * matrix without diagonal */
void
outerselfsubtr(double* x, int n, double *out) {
  int k=0;

  /* do not include diagonal */
  for (int j=0; j < n-1; j++)
    for (int i=j+1; i < n; i++)
      out[k++] = x[i] - x[j];
}

/* calculate kernel cumulative distribution function values
 * on the nonzero entries (only) from the input sparse matrix,
 * which should be provided in both, the compressed sparse column
 * (CSC) and the compressed sparse row (CSR) formats, where the latter
 * is used for speeding up the scanning of the rows of the input matrix.
 * the returned value is a sparse (CSC) matrix.
 */
SEXP
kcdfvals_sparse_to_sparse_R(SEXP XCspR, SEXP XRspR, SEXP GausskR, SEXP verboseR) {
  SEXP kcdfRobj;
  int* XCsp_dim;
  int* XCsp_i;
  int* XCsp_p;
  double* XCsp_x;
  int* XRsp_j;
  int* XRsp_p;
  double* XRsp_x;
  int  nnz = length(GET_SLOT(XCspR, Matrix_xSym));
  Rboolean Gaussk=asLogical(GausskR);
  Rboolean verbose=asLogical(verboseR);
  int* kcdfRobj_dim;
  int* kcdfRobj_i;
  int* kcdfRobj_p;
  double* kcdfRobj_x;
  int  nr, nc;
  SEXP pb = R_NilValue;
  int  nunprotect=0;

  PROTECT(XCspR); nunprotect++;
  PROTECT(XRspR); nunprotect++;

  XCsp_dim = INTEGER(GET_SLOT(XCspR, Matrix_DimSym));
  nr = XCsp_dim[0]; /* number of rows */
  nc = XCsp_dim[1]; /* number of columns */
  XCsp_i = INTEGER(GET_SLOT(XCspR, Matrix_iSym));
  XCsp_p = INTEGER(GET_SLOT(XCspR, Matrix_pSym));
  XCsp_x = REAL(GET_SLOT(XCspR, Matrix_xSym));

  XRsp_j = INTEGER(GET_SLOT(XRspR, Matrix_jSym));
  XRsp_p = INTEGER(GET_SLOT(XRspR, Matrix_pSym));
  XRsp_x = REAL(GET_SLOT(XRspR, Matrix_xSym));

  /* create a new dgCMatrix object (CSC) to store the result,
   * copying the i and p slots from the input CsparseMatrix,
   * and allocating memory for the x slot */
  kcdfRobj = PROTECT(NEW_OBJECT(MAKE_CLASS("dgCMatrix"))); nunprotect++;
  kcdfRobj_dim = INTEGER(ALLOC_SLOT(kcdfRobj, Matrix_DimSym, INTSXP, 2));
  kcdfRobj_dim[0] = nr;
  kcdfRobj_dim[1] = nc;
  kcdfRobj_i = INTEGER(ALLOC_SLOT(kcdfRobj, Matrix_iSym, INTSXP, nnz));
  kcdfRobj_p = INTEGER(ALLOC_SLOT(kcdfRobj, Matrix_pSym, INTSXP, nc+1));
  kcdfRobj_x = REAL(ALLOC_SLOT(kcdfRobj, Matrix_xSym, REALSXP, nnz));
  Memcpy(kcdfRobj_i, XCsp_i, (size_t) nnz);
  Memcpy(kcdfRobj_p, XCsp_p, (size_t) (nc+1));
  Memcpy(kcdfRobj_x, XCsp_x, (size_t) nnz);

  if (verbose) {
    pb = PROTECT(cli_progress_bar(nr, NULL));
    cli_progress_set_name(pb, "Estimating ECDFs");
    nunprotect++;
  }

  for (int i=0; i < nr; i++) {
    int nv;

    if (verbose) { /* show progress */
      if (i % 100 == 0 && CLI_SHOULD_TICK)
        cli_progress_set(pb, i);
    }

    /* number of nonzero values in the i-th row */
    nv = XRsp_p[i+1]-XRsp_p[i];
    if (nv > 0) {
      double* x = XRsp_x + XRsp_p[i];
      double* r = R_Calloc(nv, double);
      /*
      double  h;

      h = Gaussk ? sd(x, nv) : 0.5;
      ossx = R_Calloc(nv, double);
      outerselfsubtr(x, nv, ossx);

        kx = Gaussk ? pnorm(ossx[j], 0.0, 1.0, TRUE, FALSE) :
                      ppois(ossx[j], ossx[j]+h, TRUE, FALSE);
      */
      row_d_nologodds(x, x, r, nv, nv, !Gaussk);

      /* set kcdf values on the corresponding positions
       * of the output CSC matrix */
      for (int j=XRsp_p[i]; j < XRsp_p[i+1]; j++) {
        int k = j - XRsp_p[i];          /* index value at kcdf */
        int col = XRsp_j[j];            /* zero-based col index */
        int nnz2col = XCsp_p[col];      /* nonzero values to col */
        int s;                          /* shift to row in the col */
        int XCsp_idx;                   /* index value at CSC */
     
        s = 0;
        while (XCsp_i[nnz2col+s] != i && nnz2col+s < XCsp_p[col+1])
          s++;
        XCsp_idx = nnz2col + s;
        kcdfRobj_x[XCsp_idx] = r[k];
      }

      R_Free(r);
    }

  }

  if (verbose)
    cli_progress_done(pb);

  UNPROTECT(nunprotect++); /* XCspR XRspR kcdfRobj pb */

  return(kcdfRobj);
}



/* calculate kernel cumulative distribution function values
 * on the zero and nonzero entries from the input sparse matrix,
 * which should be provided in both, the compressed sparse column
 * (CSC) and the compressed sparse row (CSR) formats, where the latter
 * is used for speeding up the scanning of the rows of the input matrix.
 * the returned value is a dense matrix.
 */
SEXP
kcdfvals_sparse_to_dense_R(SEXP XCspR, SEXP XRspR, SEXP GausskR, SEXP verboseR) {
  SEXP kcdfRobj;
  double* kcdf_vals;
  int* XCsp_dim;
  int* XRsp_j;
  int* XRsp_p;
  double* XRsp_x;
  Rboolean Gaussk=asLogical(GausskR);
  Rboolean verbose=asLogical(verboseR);
  int  nr, nc;
  SEXP pb = R_NilValue;
  int  nunprotect=0;

  PROTECT(XCspR); nunprotect++;
  PROTECT(XRspR); nunprotect++;

  XCsp_dim = INTEGER(GET_SLOT(XCspR, Matrix_DimSym));
  nr = XCsp_dim[0]; /* number of rows */
  nc = XCsp_dim[1]; /* number of columns */

  XRsp_j = INTEGER(GET_SLOT(XRspR, Matrix_jSym));
  XRsp_p = INTEGER(GET_SLOT(XRspR, Matrix_pSym));
  XRsp_x = REAL(GET_SLOT(XRspR, Matrix_xSym));

  /* create a new dense matrix object to store the result,
   * if nr * nc > INT_MAX and LONG_VECTOR_SUPPORT is not
   * available, the function allocMatrix() will prompt an error */
  kcdfRobj  = PROTECT(allocMatrix(REALSXP, nr, nc)); nunprotect++;
  kcdf_vals = REAL(kcdfRobj);

  if (verbose) {
    pb = PROTECT(cli_progress_bar(nr, NULL));
    cli_progress_set_name(pb, "Estimating ECDFs");
    nunprotect++;
  }

  for (int i=0; i < nr; i++) {
    double*       x = R_Calloc(nc, double); /* assuming zeroes are set */
    double*       r = R_Calloc(nc, double); /* assuming zeroes are set */

    if (verbose) { /* show progress */
      if (i % 100 == 0 && CLI_SHOULD_TICK)
        cli_progress_set(pb, i);
    }

    /* convert sparse row into a dense vector */
    for (int j=XRsp_p[i]; j < XRsp_p[i+1]; j++)
      x[XRsp_j[j]] = XRsp_x[j];

    row_d_nologodds(x, x, r, nc, nc, !Gaussk);

    /* set kcdf values on the corresponding positions
     * of the output dense matrix */
    for (int j=0; j < nc; j++) {
#ifdef LONG_VECTOR_SUPPORT
      R_xlen_t idx = nr * j + i;
#else
      int idx = nr * j + i;
#endif

      kcdf_vals[idx] = r[j];
    }

    R_Free(r);
    R_Free(x);
  }

  if (verbose)
    cli_progress_done(pb);

  UNPROTECT(nunprotect); /* XCspR XRspR kcdfRobj pb */

  return(kcdfRobj);
}
