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
            Matrix_pSym,
            SVT_SparseArray_typeSym,
            SVT_SparseArray_dimNamesSym,
            SVT_SparseArray_dimSym,
            SVT_SparseArray_svtSym;

double
sd(double* x, int n);

void
outerselfsubtr(double* x, int n, double *out); 

void
row_d_nologodds(double* x, double* y, double* r, int size_density_n,
                int size_test_n, int Gaussk);

int
fetch_row_nzvals(SEXP svtR, int i, int itypevals, int* whimin1,
                 int* inzvals, double* dnzvals, int* nzcols);

SEXP
kcdfvals_svt_to_dense_R(SEXP XsvtR, SEXP GausskR, SEXP verboseR);

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
      row_d_nologodds(x, x, r, nv, nv, Gaussk);

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

    row_d_nologodds(x, x, r, nc, nc, Gaussk);

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

/* calculate kernel cumulative distribution function values
 * on the zero and nonzero entries from the input sparse matrix,
 * which should be provided as an SVT_SparseMatrix object.
 * the returned value is a dense matrix.
 */
SEXP
kcdfvals_svt_to_dense_R(SEXP XsvtR, SEXP GausskR, SEXP verboseR) {
  SEXP        kcdfRobj;
  double*     kcdf_vals;
  SEXP        Xsvt_dimR;
  int*        Xsvt_dim;
  SEXP        Xsvt_SVT;
  const char* Xsvt_type;
  Rboolean    Gaussk=asLogical(GausskR);
  Rboolean    verbose=asLogical(verboseR);
  int         itypevals;
  int*        nzcols;
  int*        whimin1;
  int*        inzvals=NULL;
  double*     dnzvals=NULL;
  int         nr, nc;
  SEXP        pb=R_NilValue;
  int         nunprotect=0;

  PROTECT(XsvtR); nunprotect++;

  Xsvt_dimR = GET_SLOT(XsvtR, SVT_SparseArray_dimSym);
  if (length(Xsvt_dimR) > 2)
    error("the input SVT_SparseMatrix object can only have two dimensions and has %d", length(Xsvt_dimR));
  Xsvt_dim = INTEGER(Xsvt_dimR);
  nr = Xsvt_dim[0]; /* number of rows */
  nc = Xsvt_dim[1]; /* number of columns */

  Xsvt_type = CHAR(STRING_ELT(getAttrib(XsvtR, SVT_SparseArray_typeSym), 0));
  Xsvt_SVT = GET_SLOT(XsvtR, SVT_SparseArray_svtSym);

  nzcols = R_Calloc(nc, int);
  whimin1 = R_Calloc(nc, int);
  itypevals = 0;
  if (!strcmp(Xsvt_type, "integer")) {
    itypevals = 1;
    inzvals = R_Calloc(nc, int);
  } else
    dnzvals = R_Calloc(nc, double);

  /* create a new dense matrix object to store the result,
   * if nr * nc > INT_MAX and LONG_VECTOR_SUPPORT is not
   * available, the function allocMatrix() will prompt an error */
  kcdfRobj = PROTECT(allocMatrix(REALSXP, nr, nc)); nunprotect++;
  kcdf_vals = REAL(kcdfRobj);

  if (verbose) {
    pb = PROTECT(cli_progress_bar(nr, NULL));
    cli_progress_set_name(pb, "Estimating ECDFs");
    nunprotect++;
  }

  for (int i=0; i < nr; i++) {
    int           nv;
    double*       x=R_Calloc(nc, double); /* assuming zeroes are set */
    double*       r=R_Calloc(nc, double); /* assuming zeroes are set */

    if (verbose) { /* show progress */
      if (i % 100 == 0 && CLI_SHOULD_TICK)
        cli_progress_set(pb, i);
    }

    /* fetch nonzero values in the i-th row */
    nv = fetch_row_nzvals(Xsvt_SVT, i, itypevals, whimin1,
                          inzvals, dnzvals, nzcols);

    /* convert sparse row into a dense vector */
    for (int j=0; j < nv; j++)
      x[nzcols[j]] = itypevals ? ((double) inzvals[j]) : dnzvals[j];

    row_d_nologodds(x, x, r, nc, nc, Gaussk);

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

  R_Free(nzcols);
  R_Free(whimin1);
  if (itypevals)
    R_Free(inzvals);
  else
    R_Free(dnzvals);

  if (verbose)
    cli_progress_done(pb);

  UNPROTECT(nunprotect); /* XsvtR kcdfRobj pb */

  return(kcdfRobj);
}

/* calculate kernel cumulative distribution function values
 * on the nonzero entries (only) from the input sparse matrix,
 * which should be provided as a SVT_SparseArray object.
 * the returned value is a SVT_SparseArray object.
 */
SEXP
kcdfvals_svt_to_svt_R(SEXP XsvtR, SEXP GausskR, SEXP verboseR) {
  SEXP        Xsvt_dimR;
  int*        Xsvt_dim;
  SEXP        kcdfRobj;
  SEXP        Xsvt_SVT;
  const char* Xsvt_type;
  Rboolean    Gaussk=asLogical(GausskR);
  Rboolean    verbose=asLogical(verboseR);
  int*        kcdfRobj_dim;
  int*        nnzcols; /* counter of nonzero values as they get filled up */
  SEXP        kcdfRobj_SVT;
  int         itypevals;
  int*        nzcols; /* 0-based index of the columns with nonzero values */
  int*        whimin1;
  int*        inzvals = NULL;
  double*     dnzvals = NULL;
  int         nr, nc;
  SEXP        pb=R_NilValue;
  int         nunprotect=0;

  PROTECT(XsvtR); nunprotect++;

  Xsvt_dimR = GET_SLOT(XsvtR, SVT_SparseArray_dimSym);
  if (length(Xsvt_dimR) > 2)
    error("the input SVT_SparseMatrix object can only have two dimensions and has %d", length(Xsvt_dimR));
  Xsvt_dim = INTEGER(Xsvt_dimR);
  nr = Xsvt_dim[0]; /* number of rows */
  nc = Xsvt_dim[1]; /* number of columns */
  Xsvt_type = CHAR(STRING_ELT(getAttrib(XsvtR, SVT_SparseArray_typeSym), 0));
  Xsvt_SVT = GET_SLOT(XsvtR, SVT_SparseArray_svtSym);

  itypevals = 0;
  if (!strcmp(Xsvt_type, "integer")) {
    itypevals = 1;
    inzvals = R_Calloc(nc, int);
  } else
    dnzvals = R_Calloc(nc, double);

  /* create a new SVT_SparseMatrix object to store the result */
  kcdfRobj = PROTECT(NEW_OBJECT(MAKE_CLASS("SVT_SparseMatrix"))); nunprotect++;
  kcdfRobj_dim = INTEGER(ALLOC_SLOT(kcdfRobj, SVT_SparseArray_dimSym, INTSXP, 2));
  kcdfRobj_dim[0] = nr;
  kcdfRobj_dim[1] = nc;
  SET_STRING_ELT(ALLOC_SLOT(kcdfRobj, SVT_SparseArray_typeSym, STRSXP, 1), 0,
                 mkChar("double"));
  SET_SLOT(kcdfRobj, SVT_SparseArray_svtSym, duplicate(Xsvt_SVT));
  kcdfRobj_SVT = GET_SLOT(kcdfRobj, SVT_SparseArray_svtSym);
  if (itypevals) { /* if input is integer then replace integer values by double */
    for (int i=0; i < nc; i++) {
      if (VECTOR_ELT(kcdfRobj_SVT, i) != R_NilValue)
        SET_VECTOR_ELT(VECTOR_ELT(kcdfRobj_SVT, i), 0,
                       coerceVector(VECTOR_ELT(VECTOR_ELT(kcdfRobj_SVT, i), 0),
                                    REALSXP));
    }
  }

  nnzcols = R_Calloc(nc, int); /* assuming values are initialized to 0 */

  nzcols = R_Calloc(nc, int);
  whimin1 = R_Calloc(nc, int);
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

    /* fetch nonzero values in the i-th row */
    nv = fetch_row_nzvals(Xsvt_SVT, i, itypevals, whimin1,
                          inzvals, dnzvals, nzcols);

    if (nv > 0) {
      double* x = dnzvals;
      double* r = R_Calloc(nv, double);

      if (itypevals) {
        x = R_Calloc(nv, double);
        for (int j=0; j < nv; j++)
          x[j] = (double) inzvals[j];
      }
      row_d_nologodds(x, x, r, nv, nv, Gaussk);

      /* set kcdf values on the corresponding positions
       * of the output SVT_SparseArray object */
      for (int j=0; j < nv; j++) {
        int     col = nzcols[j];            /* zero-based col index */
        double* y = REAL(VECTOR_ELT(VECTOR_ELT(kcdfRobj_SVT, col), 0));

        y[nnzcols[col]] = r[j];
        nnzcols[col]++;
      }

      if (itypevals)
        R_Free(x);

      R_Free(r);
    }

  }

  R_Free(whimin1);
  R_Free(nzcols);
  R_Free(nnzcols);
  if (itypevals)
    R_Free(inzvals);
  else
    R_Free(dnzvals);

  if (verbose)
    cli_progress_done(pb);

  UNPROTECT(nunprotect++); /* XsvtR kcdfRobj pb */

  return(kcdfRobj);
}
