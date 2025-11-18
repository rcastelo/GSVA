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

#ifdef LONG_VECTOR_SUPPORT
R_xlen_t
nzcount_intCSCp_SVT(SEXP SVT, int* CSCp);
#else
int
nzcount_intCSCp_SVT(SEXP SVT, int* CSCp);
#endif

SEXP
fetch_row_nzvals_R(SEXP svtR, SEXP iR);

int
fetch_row_nzvals(SEXP svtR, int i, int itypevals, int* inzvals,
                 double* dnzvals, int* nzcols);

SEXP
ecdfvals_svt_to_dense_R(SEXP XsvtR, SEXP verboseR);

SEXP
ecdfvals_svt_to_sparse_R(SEXP XsvtR, SEXP verboseR);

SEXP
ecdfvals_svt_to_svt_R(SEXP XsvtR, SEXP verboseR);

SEXP
ecdfvals_sparse_to_sparse_R(SEXP XCspR, SEXP XRspR, SEXP verboseR);

SEXP
ecdfvals_sparse_to_dense_R(SEXP XCspR, SEXP XRspR, SEXP verboseR);

SEXP
ecdfvals_dense_to_dense_R(SEXP XR, SEXP verboseR);

SEXP
ecdfvals_dense_to_dense_nas_R(SEXP XR, SEXP verboseR);

SEXP
match_int(SEXP x, SEXP table);

int
pending_interrupt(void);

SEXP
allocMatrix(SEXPTYPE mode, int nrow, int ncol);

/* numerical (double) comparison for qsort() - not needed with R_qsort() */
int
dbl_cmp(const void* a, const void* b) {
  const double *ad = (double*) a;
  const double *bd = (double*) b;

  if (*ad > *bd)
    return 1;
  if (*ad < *bd)
    return -1;

  return 0;
}

/* this is a wrapper to call fetch_row_nzvals() directly from R,
 * mostly for testing purposes (i.e., it is not exported) */
SEXP
fetch_row_nzvals_R(SEXP XsvtR, SEXP iR) {
  int     nr, nc;
  int*    Xsvt_dim;
  const char*   Xsvt_type;
  SEXP    Xsvt_SVT;
  int     i = INTEGER(iR)[0];
  int     itypevals;
  int*    inzvals = NULL;
  double* dnzvals = NULL;
  int*    nzcols;
  int     nnzvals;
  SEXP    ansR, nzcR, nzvR;

  Xsvt_dim = INTEGER(GET_SLOT(XsvtR, SVT_SparseArray_dimSym));
  nr = Xsvt_dim[0]; /* number of rows */
  nc = Xsvt_dim[1]; /* number of columns */
  Xsvt_type = CHAR(STRING_ELT(getAttrib(XsvtR, SVT_SparseArray_typeSym), 0));

  if (i < 1 || i > nr)
    error("i (%d) should be a number between 1 and the number of rows in the input SVT_SparseArray matrix (%d)", i, nr);

  itypevals = 0;
  if (!strcmp(Xsvt_type, "integer"))
    itypevals = 1;

  Xsvt_SVT = GET_SLOT(XsvtR, SVT_SparseArray_svtSym);

  PROTECT(ansR = allocVector(VECSXP, 2));
  PROTECT(nzcR = allocVector(INTSXP, nc));

  if (itypevals) {
    PROTECT(nzvR = allocVector(INTSXP, nc));
    inzvals = INTEGER(nzvR);
  } else {
    PROTECT(nzvR = allocVector(REALSXP, nc));
    dnzvals = REAL(nzvR);
  }
  nzcols = INTEGER(nzcR);

  nnzvals = fetch_row_nzvals(Xsvt_SVT, i-1, itypevals,  /* input  */
                             inzvals, dnzvals, nzcols); /* output */

  /* following advice in https://stat.ethz.ch/pipermail/r-package-devel/2025q1/011366.html */
  nzcR = Rf_lengthgets(nzcR, nnzvals);
  nzvR = Rf_lengthgets(nzvR, nnzvals);
  
  SET_VECTOR_ELT(ansR, 0, nzvR);
  SET_VECTOR_ELT(ansR, 1, nzcR);

  UNPROTECT(3); /* ansR nzvR nzcR */

  return(ansR);
}

/* fetch the non-zero values from given 0-based (i) row of an input
 * SVT_SparseArray slot @SVT list object (svtR). Depending on whether
 * the input values * are integer or double (itypevals), the non-zero
 * values are returned by reference into inzvals (integer) or dnzvals
 * (double). the columns of svtR having non-zero values for the row i
 * are returned by reference * in nzcols, and the function returns the
 * value of nnzvals containing the number of non-zero values in row i
 */
int
fetch_row_nzvals(SEXP svtR, int i, int itypevals, int* inzvals,
                 double* dnzvals, int* nzcols) {
  int     nc=length(svtR);
  int     nnzvals;

  nnzvals = 0;
  if (length(svtR) == 0) /* when matrix is empty, SVT is NULL */
    return(nnzvals);

  for (int j=0; j < nc; j++) {
    SEXP svtLeaf = VECTOR_ELT(svtR, j);

    if (svtLeaf != R_NilValue) {
      SEXP    valsR = VECTOR_ELT(svtLeaf, 0);
      SEXP    offsetsR = VECTOR_ELT(svtLeaf, 1);
      int     noffsets = length(offsetsR);
      int     nvals = length(valsR);
      int*    offsets = INTEGER(offsetsR);
      int*    ivals;
      double* dvals;
      int     k;

      if (itypevals)
        ivals = INTEGER(valsR);
      else
        dvals = REAL(valsR);

      k=0;
      while (k < noffsets && offsets[k] < i)
        k++;

      if (k < noffsets) {
        if (offsets[k] == i) {
          if (itypevals)
            inzvals[nnzvals] = nvals > 0 ? ivals[k] : 1;
          else
            dnzvals[nnzvals] = nvals > 0 ? dvals[k] : 1;
          if (nzcols != NULL)
            nzcols[nnzvals] = j;
          nnzvals++;
        }
      }
    }
  }

  return(nnzvals);
}

/* calculate empirical cumulative distribution function values
 * on the zero and nonzero entries from the input sparse matrix,
 * provided as an SVT_SparseMatrix object.
 * the returned value is a dense matrix.
 */
SEXP
ecdfvals_svt_to_dense_R(SEXP XsvtR, SEXP verboseR) {
  SEXP        ecdfRobj;
  double*     ecdf_vals;
  SEXP        Xsvt_dimR;
  int*        Xsvt_dim;
  SEXP        Xsvt_SVT;
  const char* Xsvt_type;
  Rboolean    verbose=asLogical(verboseR);
  int         itypevals;
  int*        nzcols;
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

  nzcols = R_Calloc(nc, int);
  itypevals = 0;
  if (!strcmp(Xsvt_type, "integer")) {
    itypevals = 1;
    inzvals = R_Calloc(nc, int);
  } else
    dnzvals = R_Calloc(nc, double);

  /* create a new dense matrix object to store the result,
   * if nr * nc > INT_MAX and LONG_VECTOR_SUPPORT is not
   * available, the function allocMatrix() will prompt an error */
  ecdfRobj = PROTECT(allocMatrix(REALSXP, nr, nc)); nunprotect++;

  if (verbose) {
    pb = PROTECT(cli_progress_bar(nr, NULL));
    cli_progress_set_name(pb, "Estimating ECDFs");
    nunprotect++;
  }

  for (int i=0; i < nr; i++) {
    SEXP          uniqvR;
    int           nv, nuniqv;
    SEXP          xR;
    double*       x;
    double*       uniqv;
    double*       ecdfuniqv;
    int           sum;
    double*       e1_p;
    const double* e2_p;
    int*          mt;
    int*          tab;
    Rboolean      zeroes=FALSE;
    int           whz, icz;

    if (verbose) { /* show progress */
      if (i % 100 == 0 && CLI_SHOULD_TICK)
        cli_progress_set(pb, i);
    }

    /* fetch nonzero values in the i-th row */
    nv = fetch_row_nzvals(Xsvt_SVT, i, itypevals, inzvals, dnzvals, nzcols);

    if (nv < nc) { /* if there is at least one zero in the row */
      nv++;        /* add that zero as an extra possible value */
      zeroes=TRUE;
    }

    /* remove consecutive repeated elements */
    /* consider adding LONG_VECTOR_SUPPORT */
    PROTECT(uniqvR = allocVector(REALSXP, nv));
    uniqv = REAL(uniqvR);
    PROTECT(xR = allocVector(REALSXP, zeroes ? nv-1 : nv));
    x = REAL(xR);

    if (zeroes) {   /* if there is at least one zero in the row */
      uniqv[0] = 0; /* add that zero as an extra possible value */
      for (int j=1; j < nv; j++) {
        x[j-1] = itypevals ? ((double) inzvals[j-1]) : dnzvals[j-1];
        uniqv[j] = x[j-1];
      }
    } else {
      for (int j=0; j < nv; j++) {
        uniqv[j] = x[j] = itypevals ? ((double) inzvals[j]) : dnzvals[j];
      }
    }

    R_qsort(uniqv, (size_t) 1, (size_t) nv);
    e1_p = uniqv;
    e2_p = e1_p + 1;
    nuniqv = 0;
    if (nv > 0)
      nuniqv = 1;
    for (int j=0; j < nv-1; j++) { /* -1 for e2_p = e1_p + 1 */
      if (*e2_p != *e1_p) {
        *(++e1_p) = *e2_p;
        nuniqv++;
      }
      e2_p++;
    }

    /* match original values to sorted unique values */
    /* consider adding LONG_VECTOR_SUPPORT */
    mt = INTEGER(match_int(xR, uniqvR)); /* 1-based! */

    /* tabulate matches */
    /* consider adding LONG_VECTOR_SUPPORT */
    tab = R_Calloc(nuniqv, int); /* assuming zeroes are set */
    for (int j=0; j < length(xR); j++) {
      if (mt[j] > 0 && mt[j] <= nuniqv)
        tab[mt[j] - 1]++;
    }

    whz = -1;
    if (zeroes) { /* if there is at least one zero in the row */
      int j = 0;
      while (j < nuniqv && tab[j] != 0)
        j++;
      if (j < nuniqv && tab[j] == 0) {  /* add the number of zeroes */
        tab[j] = nc - nv + 1; /* +1 b/c one zero is in nv */
        whz = j;
      }
    }

    /* cumulative sum to calculate ecdf values */
    /* consider adding LONG_VECTOR_SUPPORT */
    ecdfuniqv = R_Calloc(nuniqv, double); /* assuming zeroes are set */
    sum = 0;
    for (int j=0; j < nuniqv; j++) {
      sum = sum + tab[j];
      ecdfuniqv[j] = ((double) sum) / ((double) nc);
    }

    /* set ecdf values on the corresponding positions
     * of the output dense matrix */
    ecdf_vals = REAL(ecdfRobj);
    icz = 0; /* zero-based index of the columns at zeroes */
    for (int j=0; j < length(xR); j++) {
      int col = nzcols[j];            /* zero-based col index */
#ifdef LONG_VECTOR_SUPPORT
      R_xlen_t idx = nr * col + i;
#else
      int idx = nr * col + i;
#endif
      while (icz < col) {   /* fill up the zero columns */
#ifdef LONG_VECTOR_SUPPORT
        R_xlen_t idxz = nr * icz + i;
#else
        int idxz = nr * icz + i;
#endif
        ecdf_vals[idxz] = ecdfuniqv[whz];
        icz++;
      }
      icz = col+1;
      ecdf_vals[idx] = ecdfuniqv[mt[j]-1];
    }

    for (int j=icz; j < nc; j++) { /* fill up remaining zero columns */
#ifdef LONG_VECTOR_SUPPORT
        R_xlen_t idxz = nr * j + i;
#else
        int idxz = nr * j + i;
#endif
        ecdf_vals[idxz] = ecdfuniqv[whz];
    }

    R_Free(ecdfuniqv);
    R_Free(tab);

    UNPROTECT(2); /* xR uniqvR */
  }

  R_Free(nzcols);
  if (itypevals)
    R_Free(inzvals);
  else
    R_Free(dnzvals);

  if (verbose)
    cli_progress_done(pb);

  UNPROTECT(nunprotect); /* XsvtR ecdfRobj pb */

  return(ecdfRobj);
}

/* counts the total number of nonzero values in a 2-dimensional SVT_SparseArray
 * object, as well as its running sum returned by reference in CSCp */
#ifdef LONG_VECTOR_SUPPORT
R_xlen_t
nzcount_intCSCp_SVT(SEXP SVT, int* CSCp) {
   R_xlen_t nzcount = 0;
#else
int
nzcount_intCSCp_SVT(SEXP SVT, int* CSCp) {
  int nzcount = 0;
#endif
  int SVT_len = length(SVT);
  CSCp[0] = 0;
  for (int i=0; i < SVT_len; i++) {
    SEXP subSVT = VECTOR_ELT(SVT, i);

#ifdef LONG_VECTOR_SUPPORT
    nzcount = nzcount + XLENGTH(VECTOR_ELT(subSVT, 1));
#else
    nzcount = nzcount + length(VECTOR_ELT(subSVT, 1));
#endif
    CSCp[i+1] = (int) nzcount;
  }

  return(nzcount);
}

/* calculate empirical cumulative distribution function values
 * on the nonzero entries (only) from the input sparse matrix,
 * which should be provided as a SVT_SparseArray object.
 * the returned value is a sparse (CSC) dgCMatrix object.
 */
SEXP
ecdfvals_svt_to_sparse_R(SEXP XsvtR, SEXP verboseR) {
  SEXP        Xsvt_dimR;
  int*        Xsvt_dim;
#ifdef LONG_VECTOR_SUPPORT
  R_xlen_t    nnz;
#else
  int         nnz;
#endif
  SEXP        ecdfRobj;
  SEXP        Xsvt_SVT;
  const char* Xsvt_type;
  Rboolean    verbose=asLogical(verboseR);
  int*        ecdfRobj_dim;
  int*        ecdfRobj_i;
  int*        nnzcols; /* counter of nonzero values as they get filled up */
  int*        ecdfRobj_p;
  double*     ecdfRobj_x;
  int         itypevals;
  int*        nzcols; /* 0-based index of the columns with nonzero values */
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

  Xsvt_SVT = GET_SLOT(XsvtR, SVT_SparseArray_svtSym);

  /* create a new dgCMatrix object (CSC) to store the result */
  ecdfRobj = PROTECT(NEW_OBJECT(MAKE_CLASS("dgCMatrix"))); nunprotect++;
  ecdfRobj_dim = INTEGER(ALLOC_SLOT(ecdfRobj, Matrix_DimSym, INTSXP, 2));
  ecdfRobj_dim[0] = nr;
  ecdfRobj_dim[1] = nc;
  ecdfRobj_p = INTEGER(ALLOC_SLOT(ecdfRobj, Matrix_pSym, INTSXP, nc+1));

  nnz = nzcount_intCSCp_SVT(Xsvt_SVT, ecdfRobj_p);

  if (nnz > INT_MAX) {
    UNPROTECT(nunprotect);

    error("input SVT_SparseArray matrix has more non-zero values than a dgCMatrix object can store");
  }

  Xsvt_type = CHAR(STRING_ELT(getAttrib(XsvtR, SVT_SparseArray_typeSym), 0));
  Xsvt_SVT = GET_SLOT(XsvtR, SVT_SparseArray_svtSym);

  nnzcols = R_Calloc(nc, int); /* assuming values are initialized to 0 */
  nzcols = R_Calloc(nc, int);
  itypevals = 0;
  if (!strcmp(Xsvt_type, "integer")) {
    itypevals = 1;
    inzvals = R_Calloc(nc, int);
  } else
    dnzvals = R_Calloc(nc, double);

  /* create a new dgCMatrix object (CSC) to store the result (cont'ed) */
  ecdfRobj_i = INTEGER(ALLOC_SLOT(ecdfRobj, Matrix_iSym, INTSXP, nnz));
  ecdfRobj_x = REAL(ALLOC_SLOT(ecdfRobj, Matrix_xSym, REALSXP, nnz));

  if (verbose) {
    pb = PROTECT(cli_progress_bar(nr, NULL));
    cli_progress_set_name(pb, "Estimating ECDFs");
    nunprotect++;
  }

  for (int i=0; i < nr; i++) {
    SEXP          xR, uniqvR;
    int           nv, nuniqv;
    double*       x;
    double*       uniqv;
    double*       ecdfuniqv;
    int           sum;
    double*       e1_p;
    const double* e2_p;
    int*          mt;
    int*          tab;

    if (verbose) { /* show progress */
      if (i % 100 == 0 && CLI_SHOULD_TICK)
        cli_progress_set(pb, i);
    }

    /* fetch nonzero values in the i-th row */
    nv = fetch_row_nzvals(Xsvt_SVT, i, itypevals, inzvals, dnzvals, nzcols);

    /* remove consecutive repeated elements */
    /* consider adding LONG_VECTOR_SUPPORT */
    PROTECT(uniqvR = allocVector(REALSXP, nv));
    uniqv = REAL(uniqvR);
    PROTECT(xR = allocVector(REALSXP, nv));
    x = REAL(xR);
    for (int j=0; j < nv; j++) {
      uniqv[j] = x[j] = itypevals ? ((double) inzvals[j]) : dnzvals[j];
    }

    R_qsort(uniqv, (size_t) 1, (size_t) nv);
    e1_p = uniqv;
    e2_p = e1_p + 1;
    nuniqv = 0;
    if (nv > 0)
      nuniqv = 1;
    for (int j=0; j < nv-1; j++) { /* -1 for e2_p = e1_p + 1 */
      if (*e2_p != *e1_p) {
        *(++e1_p) = *e2_p;
        nuniqv++;
      }
      e2_p++;
    }

    /* match original values to sorted unique values */
    /* consider adding LONG_VECTOR_SUPPORT */
    mt = INTEGER(match_int(xR, uniqvR)); /* 1-based! */

    /* tabulate matches */
    /* consider adding LONG_VECTOR_SUPPORT */
    tab = R_Calloc(nuniqv, int); /* assuming zeroes are set */
    for (int j=0; j < length(xR); j++) {
      if (mt[j] > 0 && mt[j] <= nuniqv)
        tab[mt[j] - 1]++;
    }

    /* cumulative sum to calculate ecdf values */
    /* consider adding LONG_VECTOR_SUPPORT */
    ecdfuniqv = R_Calloc(nuniqv, double); /* assuming zeroes are set */
    sum = 0;
    for (int j=0; j < nuniqv; j++) {
      sum = sum + tab[j];
      ecdfuniqv[j] = ((double) sum) / ((double) nv);
    }

    /* set ecdf values on the corresponding positions
     * of the output CSC matrix */
    for (int j=0; j < nv; j++) {
      int col = nzcols[j];                      /* zero-based col index */
      int idx = ecdfRobj_p[col] + nnzcols[col];
    
      ecdfRobj_i[idx] = i;
      ecdfRobj_x[idx] = ecdfuniqv[mt[j]-1];
      nnzcols[col]++;
    }

    R_Free(ecdfuniqv);
    R_Free(tab);

    UNPROTECT(2); /* xR uniqvR */
  }

  R_Free(nzcols);
  R_Free(nnzcols);
  if (itypevals)
    R_Free(inzvals);
  else
    R_Free(dnzvals);

  if (verbose)
    cli_progress_done(pb);

  UNPROTECT(nunprotect); /* XsvtR ecdfRobj pb */

  return(ecdfRobj);
}


/* calculate empirical cumulative distribution function values
 * on the nonzero entries (only) from the input sparse matrix,
 * which should be provided as a SVT_SparseArray object.
 * the returned value is a SVT_SparseArray object.
 */
SEXP
ecdfvals_svt_to_svt_R(SEXP XsvtR, SEXP verboseR) {
  SEXP        Xsvt_dimR;
  int*        Xsvt_dim;
  SEXP        ecdfRobj;
  SEXP        Xsvt_SVT;
  const char* Xsvt_type;
  Rboolean    verbose=asLogical(verboseR);
  int*        ecdfRobj_dim;
  int*        nnzcols; /* counter of nonzero values as they get filled up */
  SEXP        ecdfRobj_SVT;
  int         itypevals;
  int*        nzcols; /* 0-based index of the columns with nonzero values */
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

  itypevals = 0;
  if (!strcmp(Xsvt_type, "integer")) {
    itypevals = 1;
    inzvals = R_Calloc(nc, int);
  } else
    dnzvals = R_Calloc(nc, double);

  /* create a new SVT_SparseMatrix object to store the result */
  ecdfRobj = PROTECT(NEW_OBJECT(MAKE_CLASS("SVT_SparseMatrix"))); nunprotect++;
  ecdfRobj_dim = INTEGER(ALLOC_SLOT(ecdfRobj, SVT_SparseArray_dimSym, INTSXP, 2));
  ecdfRobj_dim[0] = nr;
  ecdfRobj_dim[1] = nc;
  SET_STRING_ELT(ALLOC_SLOT(ecdfRobj, SVT_SparseArray_typeSym, STRSXP, 1), 0,
                 mkChar("double"));
  SET_SLOT(ecdfRobj, SVT_SparseArray_svtSym, duplicate(Xsvt_SVT));
  ecdfRobj_SVT = GET_SLOT(ecdfRobj, SVT_SparseArray_svtSym);
  if (itypevals) { /* if input is integer then replace integer values by double */
    for (int i=0; i < nc; i++) {
      if (VECTOR_ELT(ecdfRobj_SVT, i) != R_NilValue)
        SET_VECTOR_ELT(VECTOR_ELT(ecdfRobj_SVT, i), 0,
                       coerceVector(VECTOR_ELT(VECTOR_ELT(ecdfRobj_SVT, i), 0),
                                    REALSXP));
    }
  }

  nnzcols = R_Calloc(nc, int); /* assuming values are initialized to 0 */
 
  nzcols = R_Calloc(nc, int);
  if (verbose) {
    pb = PROTECT(cli_progress_bar(nr, NULL));
    cli_progress_set_name(pb, "Estimating ECDFs");
    nunprotect++;
  }

  for (int i=0; i < nr; i++) {
    SEXP          xR, uniqvR;
    int           nv, nuniqv;
    double*       x;
    double*       uniqv;
    double*       ecdfuniqv;
    int           sum;
    double*       e1_p;
    const double* e2_p;
    int*          mt;
    int*          tab;

    if (verbose) { /* show progress */
      if (i % 100 == 0 && CLI_SHOULD_TICK)
        cli_progress_set(pb, i);
    }

    /* fetch nonzero values in the i-th row */
    nv = fetch_row_nzvals(Xsvt_SVT, i, itypevals, inzvals, dnzvals, nzcols);

    /* remove consecutive repeated elements */
    /* consider adding LONG_VECTOR_SUPPORT */
    PROTECT(uniqvR = allocVector(REALSXP, nv));
    uniqv = REAL(uniqvR);
    PROTECT(xR = allocVector(REALSXP, nv));
    x = REAL(xR);
    for (int j=0; j < nv; j++) {
      uniqv[j] = x[j] = itypevals ? ((double) inzvals[j]) : dnzvals[j];
    }

    R_qsort(uniqv, (size_t) 1, (size_t) nv);
    e1_p = uniqv;
    e2_p = e1_p + 1;
    nuniqv = 0;
    if (nv > 0)
      nuniqv = 1;
    for (int j=0; j < nv-1; j++) { /* -1 for e2_p = e1_p + 1 */
      if (*e2_p != *e1_p) {
        *(++e1_p) = *e2_p;
        nuniqv++;
      }
      e2_p++;
    }

    /* match original values to sorted unique values */
    /* consider adding LONG_VECTOR_SUPPORT */
    mt = INTEGER(match_int(xR, uniqvR)); /* 1-based! */

    /* tabulate matches */
    /* consider adding LONG_VECTOR_SUPPORT */
    tab = R_Calloc(nuniqv, int); /* assuming zeroes are set */
    for (int j=0; j < length(xR); j++) {
      if (mt[j] > 0 && mt[j] <= nuniqv)
        tab[mt[j] - 1]++;
    }

    /* cumulative sum to calculate ecdf values */
    /* consider adding LONG_VECTOR_SUPPORT */
    ecdfuniqv = R_Calloc(nuniqv, double); /* assuming zeroes are set */
    sum = 0;
    for (int j=0; j < nuniqv; j++) {
      sum = sum + tab[j];
      ecdfuniqv[j] = ((double) sum) / ((double) nv);
    }

    /* set ecdf values on the corresponding positions
     * of the output SVT_SparseArray object */

    for (int j=0; j < nv; j++) {
      int col = nzcols[j];                       /* zero-based col index */
      double* x = REAL(VECTOR_ELT(VECTOR_ELT(ecdfRobj_SVT, col), 0));

      x[nnzcols[col]] = ecdfuniqv[mt[j]-1];
      nnzcols[col]++;
    }

    R_Free(ecdfuniqv);
    R_Free(tab);

    UNPROTECT(2); /* xR uniqvR */
  }

  R_Free(nzcols);
  R_Free(nnzcols);
  if (itypevals)
    R_Free(inzvals);
  else
    R_Free(dnzvals);

  if (verbose)
    cli_progress_done(pb);

  UNPROTECT(nunprotect); /* XsvtR ecdfRobj pb */

  return(ecdfRobj);
}


/* calculate empirical cumulative distribution function values
 * on the nonzero entries (only) from the input sparse matrix,
 * which should be provided in both, the compressed sparse column
 * (CSC) and the compressed sparse row (CSR) formats, where the latter
 * is used for speeding up the scanning of the rows of the input matrix.
 * the returned value is a sparse (CSC) matrix.
 */
SEXP
ecdfvals_sparse_to_sparse_R(SEXP XCspR, SEXP XRspR, SEXP verboseR) {
  SEXP ecdfRobj;
  int* XCsp_dim;
  int* XCsp_i;
  int* XCsp_p;
  double* XCsp_x;
  int* XRsp_j;
  int* XRsp_p;
  double* XRsp_x;
  int  nnz = length(GET_SLOT(XCspR, Matrix_xSym));
  Rboolean verbose=asLogical(verboseR);
  int* ecdfRobj_dim;
  int* ecdfRobj_i;
  int* ecdfRobj_p;
  double* ecdfRobj_x;
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
  ecdfRobj = PROTECT(NEW_OBJECT(MAKE_CLASS("dgCMatrix"))); nunprotect++;
  ecdfRobj_dim = INTEGER(ALLOC_SLOT(ecdfRobj, Matrix_DimSym, INTSXP, 2));
  ecdfRobj_dim[0] = nr;
  ecdfRobj_dim[1] = nc;
  ecdfRobj_i = INTEGER(ALLOC_SLOT(ecdfRobj, Matrix_iSym, INTSXP, nnz));
  ecdfRobj_p = INTEGER(ALLOC_SLOT(ecdfRobj, Matrix_pSym, INTSXP, nc+1));
  ecdfRobj_x = REAL(ALLOC_SLOT(ecdfRobj, Matrix_xSym, REALSXP, nnz));
  Memcpy(ecdfRobj_i, XCsp_i, (size_t) nnz);
  Memcpy(ecdfRobj_p, XCsp_p, (size_t) (nc+1));
  Memcpy(ecdfRobj_x, XCsp_x, (size_t) nnz);

  if (verbose) {
    pb = PROTECT(cli_progress_bar(nr, NULL));
    cli_progress_set_name(pb, "Estimating ECDFs");
    nunprotect++;
  }

  for (int i=0; i < nr; i++) {
    SEXP          xR, uniqvR;
    int           nv, nuniqv;
    double*       x;
    double*       uniqv;
    double*       ecdfuniqv;
    int           sum;
    double*       e1_p;
    const double* e2_p;
    int*          mt;
    int*          tab;

    if (verbose) { /* show progress */
      if (i % 100 == 0 && CLI_SHOULD_TICK)
        cli_progress_set(pb, i);
    }

    /* number of nonzero values in the i-th row */
    nv = XRsp_p[i+1]-XRsp_p[i];

    /* remove consecutive repeated elements */
    /* consider adding LONG_VECTOR_SUPPORT */
    PROTECT(uniqvR = allocVector(REALSXP, nv));
    PROTECT(xR = allocVector(REALSXP, nv));
    uniqv = REAL(uniqvR);
    x = REAL(xR);
    for (int j=XRsp_p[i]; j < XRsp_p[i+1]; j++) {
      int k = j - XRsp_p[i];
      uniqv[k] = x[k] = XRsp_x[j];
    }

    /* qsort(uniqv, nv, sizeof(double), dbl_cmp); */
    R_qsort(uniqv, (size_t) 1, (size_t) nv);
    e1_p = uniqv;
    e2_p = e1_p + 1;
    nuniqv = 0;
    if (nv > 0)
      nuniqv = 1;
    for (int j=0; j < nv-1; j++) { /* -1 for e2_p = e1_p + 1 */
      if (*e2_p != *e1_p) {
        *(++e1_p) = *e2_p;
        nuniqv++;
      }
      e2_p++;
    }

    /* match original values to sorted unique values */
    /* consider adding LONG_VECTOR_SUPPORT */
    mt = INTEGER(match_int(xR, uniqvR)); /* 1-based! */

    /* tabulate matches */
    /* consider adding LONG_VECTOR_SUPPORT */
    tab = R_Calloc(nuniqv, int); /* assuming zeroes are set */
    for (int j=XRsp_p[i]; j < XRsp_p[i+1]; j++) {
      int k = j - XRsp_p[i];
      if (mt[k] > 0 && mt[k] <= nuniqv)
        tab[mt[k] - 1]++;
    }

    /* cumulative sum to calculate ecdf values */
    /* consider adding LONG_VECTOR_SUPPORT */
    ecdfuniqv = R_Calloc(nuniqv, double); /* assuming zeroes are set */
    sum = 0;
    for (int j=0; j < nuniqv; j++) {
      sum = sum + tab[j];
      ecdfuniqv[j] = ((double) sum) / ((double) nv);
    }

    /* set ecdf values on the corresponding positions
     * of the output CSC matrix */
    for (int j=XRsp_p[i]; j < XRsp_p[i+1]; j++) {
      int k = j - XRsp_p[i];          /* index value at ecdf */
      int col = XRsp_j[j];            /* zero-based col index */
      int nnz2col = XCsp_p[col];      /* nonzero values to col */
      int s;                          /* shift to row in the col */
      int XCsp_idx;                   /* index value at CSC */
     
      s = 0;
      while (XCsp_i[nnz2col+s] != i && nnz2col+s < XCsp_p[col+1])
        s++;
      XCsp_idx = nnz2col + s;
      ecdfRobj_x[XCsp_idx] = ecdfuniqv[mt[k]-1];
    }

    R_Free(ecdfuniqv);
    R_Free(tab);

    UNPROTECT(2); /* xR uniqvR */
  }

  if (verbose)
    cli_progress_done(pb);

  UNPROTECT(nunprotect); /* XCspR XRspR ecdfRobj pb */

  return(ecdfRobj);
}


/* calculate empirical cumulative distribution function values
 * on the zero and nonzero entries from the input sparse matrix,
 * which should be provided in both, the compressed sparse column
 * (CSC) and the compressed sparse row (CSR) formats, where the latter
 * is used for speeding up the scanning of the rows of the input matrix.
 * the returned value is a dense matrix.
 */
SEXP
ecdfvals_sparse_to_dense_R(SEXP XCspR, SEXP XRspR, SEXP verboseR) {
  SEXP ecdfRobj;
  double* ecdf_vals;
  int* XCsp_dim;
  int* XRsp_j;
  int* XRsp_p;
  Rboolean verbose=asLogical(verboseR);
  double* XRsp_x;
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
  ecdfRobj = PROTECT(allocMatrix(REALSXP, nr, nc)); nunprotect++;

  if (verbose) {
    pb = PROTECT(cli_progress_bar(nr, NULL));
    cli_progress_set_name(pb, "Estimating ECDFs");
    nunprotect++;
  }

  for (int i=0; i < nr; i++) {
    SEXP          xR, uniqvR;
    int           nv, nuniqv;
    double*       x;
    double*       uniqv;
    double*       ecdfuniqv;
    int           sum;
    double*       e1_p;
    const double* e2_p;
    int*          mt;
    int*          tab;
    Rboolean      zeroes=FALSE;
    int           whz, icz;

    if (verbose) { /* show progress */
      if (i % 100 == 0 && CLI_SHOULD_TICK)
        cli_progress_set(pb, i);
    }

    /* number of nonzero values in the i-th row */
    nv = XRsp_p[i+1]-XRsp_p[i];
    if (nv < nc) { /* if there is at least one zero in the row */
      nv++;        /* add that zero as an extra possible value */
      zeroes=TRUE;
    }

    /* remove consecutive repeated elements */
    /* consider adding LONG_VECTOR_SUPPORT */
    PROTECT(uniqvR = allocVector(REALSXP, nv));
    PROTECT(xR = allocVector(REALSXP, zeroes ? nv-1 : nv));
    uniqv = REAL(uniqvR);
    x = REAL(xR);
    if (zeroes) {   /* if there is at least one zero in the row */
      uniqv[0] = 0; /* add that zero as an extra possible value */
      for (int j=XRsp_p[i]; j < XRsp_p[i+1]; j++) {
        int k = j - XRsp_p[i];
        uniqv[k+1] = XRsp_x[j];
        x[k] = XRsp_x[j];
      }
    } else {
      for (int j=XRsp_p[i]; j < XRsp_p[i+1]; j++) {
        int k = j - XRsp_p[i];
        uniqv[k] = x[k] = XRsp_x[j];
      }
    }

    /* qsort(uniqv, nv, sizeof(double), dbl_cmp); */
    R_qsort(uniqv, (size_t) 1, (size_t) nv);
    e1_p = uniqv;
    e2_p = e1_p + 1;
    nuniqv = 0;
    if (nv > 0)
      nuniqv = 1;
    for (int j=0; j < nv-1; j++) { /* -1 for e2_p = e1_p + 1 */
      if (*e2_p != *e1_p) {
        *(++e1_p) = *e2_p;
        nuniqv++;
      }
      e2_p++;
    }

    /* match original values to sorted unique values */
    /* consider adding LONG_VECTOR_SUPPORT */
    mt = INTEGER(match_int(xR, uniqvR)); /* 1-based! */

    /* tabulate matches */
    /* consider adding LONG_VECTOR_SUPPORT */
    tab = R_Calloc(nuniqv, int); /* assuming zeroes are set */
    for (int j=XRsp_p[i]; j < XRsp_p[i+1]; j++) {
      int k = j - XRsp_p[i];
      if (mt[k] > 0 && mt[k] <= nuniqv)
        tab[mt[k] - 1]++;
    }
    whz = -1;
    if (zeroes) { /* if there is at least one zero in the row */
      int j = 0;
      while (j < nuniqv && tab[j] != 0)
        j++;
      if (j < nuniqv && tab[j] == 0) {  /* add the number of zeroes */
        tab[j] = nc - nv + 1; /* +1 b/c one zero is in nv */
        whz = j;
      }
    }

    /* cumulative sum to calculate ecdf values */
    /* consider adding LONG_VECTOR_SUPPORT */
    ecdfuniqv = R_Calloc(nuniqv, double); /* assuming zeroes are set */
    sum = 0;
    for (int j=0; j < nuniqv; j++) {
      sum = sum + tab[j];
      ecdfuniqv[j] = ((double) sum) / ((double) nc);
    }

    /* set ecdf values on the corresponding positions
     * of the output dense matrix */
    ecdf_vals = REAL(ecdfRobj);
    icz = 0; /* zero-based index of the columns at zeroes */
    for (int j=XRsp_p[i]; j < XRsp_p[i+1]; j++) {
      int k = j - XRsp_p[i];          /* index value at ecdf */
      int col = XRsp_j[j];            /* zero-based col index */
#ifdef LONG_VECTOR_SUPPORT
      R_xlen_t idx = nr * col + i;
#else
      int idx = nr * col + i;
#endif
      while (icz < col) {   /* fill up the zero columns */
#ifdef LONG_VECTOR_SUPPORT
        R_xlen_t idxz = nr * icz + i;
#else
        int idxz = nr * icz + i;
#endif
        ecdf_vals[idxz] = ecdfuniqv[whz];
        icz++;
      }
      icz = col+1;
      ecdf_vals[idx] = ecdfuniqv[mt[k]-1];
    }
    for (int j=icz; j < nc; j++) { /* fill up remaining zero columns */
#ifdef LONG_VECTOR_SUPPORT
        R_xlen_t idxz = nr * j + i;
#else
        int idxz = nr * j + i;
#endif
        ecdf_vals[idxz] = ecdfuniqv[whz];
    }

    R_Free(ecdfuniqv);
    R_Free(tab);

    UNPROTECT(2); /* xR uniqvR */
  }

  if (verbose)
    cli_progress_done(pb);

  UNPROTECT(nunprotect); /* XCspR XRspR ecdfRobj pb */

  return(ecdfRobj);
}


/* calculate empirical cumulative distribution function values
 * on the zero and nonzero entries from the input dense matrix,
 * the returned value is a dense matrix.
 */
SEXP
ecdfvals_dense_to_dense_R(SEXP XR, SEXP verboseR) {
  double* X;
  Rboolean verbose=asLogical(verboseR);
  int     nr, nc;
  SEXP    ecdfRobj;
  double* ecdf_vals;
  SEXP pb = R_NilValue;
  int  nunprotect=0;

  PROTECT(XR); nunprotect++;

  nr = INTEGER(getAttrib(XR, R_DimSymbol))[0]; /* number of rows */
  nc = INTEGER(getAttrib(XR, R_DimSymbol))[1]; /* number of columns */
  X  = REAL(XR);

  /* create a new dense matrix object to store the result,
   * if nr * nc > INT_MAX and LONG_VECTOR_SUPPORT is not
   * available, the function allocMatrix() will prompt an error */
  ecdfRobj = PROTECT(allocMatrix(REALSXP, nr, nc)); nunprotect++;

  if (verbose) {
    pb = PROTECT(cli_progress_bar(nr, NULL));
    cli_progress_set_name(pb, "Estimating ECDFs");
    nunprotect++;
  }

  for (int i=0; i < nr; i++) {
    SEXP          xR, uniqvR;
    int           nuniqv;
    double*       x;
    double*       uniqv;
    double*       ecdfuniqv;
    int           sum;
    double*       e1_p;
    const double* e2_p;
    int*          mt;
    int*          tab;

    if (verbose) { /* show progress */
      if (i % 100 == 0 && CLI_SHOULD_TICK)
        cli_progress_set(pb, i);
    }

    /* remove consecutive repeated elements */
    PROTECT(uniqvR = allocVector(REALSXP, nc));
    PROTECT(xR = allocVector(REALSXP, nc));
    uniqv = REAL(uniqvR);
    x = REAL(xR);
    for (int j=0; j < nc; j++) {
#ifdef LONG_VECTOR_SUPPORT
      R_xlen_t idx = nr * j + i;
#else
      int idx = nr * j + i;
#endif
      uniqv[j] = x[j] = X[idx];
    }

    /* qsort(uniqv, nc, sizeof(double), dbl_cmp); */
    R_qsort(uniqv, (size_t) 1, (size_t) nc);
    e1_p = uniqv;
    e2_p = e1_p + 1;
    nuniqv = 0;
    if (nc > 0)
      nuniqv = 1;
    for (int j=0; j < nc-1; j++) { /* -1 for e2_p = e1_p + 1 */
      if (*e2_p != *e1_p) {
        *(++e1_p) = *e2_p;
        nuniqv++;
      }
      e2_p++;
    }

    /* match original values to sorted unique values */
    /* consider adding LONG_VECTOR_SUPPORT */
    mt = INTEGER(match_int(xR, uniqvR)); /* 1-based! */

    /* tabulate matches */
    /* consider adding LONG_VECTOR_SUPPORT */
    tab = R_Calloc(nuniqv, int); /* assuming zeroes are set */
    for (int j=0; j < nc; j++)
      if (mt[j] > 0 && mt[j] <= nuniqv)
        tab[mt[j] - 1]++;

    /* cumulative sum to calculate ecdf values */
    /* consider adding LONG_VECTOR_SUPPORT */
    ecdfuniqv = R_Calloc(nuniqv, double); /* assuming zeroes are set */
    sum = 0;
    for (int j=0; j < nuniqv; j++) {
      sum = sum + tab[j];
      ecdfuniqv[j] = ((double) sum) / ((double) nc);
    }

    /* set ecdf values on the corresponding positions
     * of the output dense matrix */
    ecdf_vals = REAL(ecdfRobj);

    for (int j=0; j < nc; j++) {
#ifdef LONG_VECTOR_SUPPORT
      R_xlen_t idx = nr * j + i;
#else
      int idx = nr * j + i;
#endif
      ecdf_vals[idx] = ecdfuniqv[mt[j]-1];
    }

    R_Free(ecdfuniqv);
    R_Free(tab);

    UNPROTECT(2); /* xR uniqvR */
  }

  if (verbose)
    cli_progress_done(pb);

  UNPROTECT(nunprotect); /* XR ecdfRobj pb */

  return(ecdfRobj);
}

/* calculate empirical cumulative distribution function values
 * on the zero and nonzero entries from the input dense matrix,
 * when missing (NA) values are present, the returned value is
 * a dense matrix.
 */
SEXP
ecdfvals_dense_to_dense_nas_R(SEXP XR, SEXP verboseR) {
  double* X;
  Rboolean verbose=asLogical(verboseR);
  int     nr, nc;
  SEXP    ecdfRobj;
  double* ecdf_vals;
  SEXP pb = R_NilValue;
  int  nunprotect=0;

  PROTECT(XR); nunprotect++;

  nr = INTEGER(getAttrib(XR, R_DimSymbol))[0]; /* number of rows */
  nc = INTEGER(getAttrib(XR, R_DimSymbol))[1]; /* number of columns */
  X  = REAL(XR);

  /* create a new dense matrix object to store the result,
   * if nr * nc > INT_MAX and LONG_VECTOR_SUPPORT is not
   * available, the function allocMatrix() will prompt an error */
  ecdfRobj = PROTECT(allocMatrix(REALSXP, nr, nc)); nunprotect++;

  if (verbose) {
    pb = PROTECT(cli_progress_bar(nr, NULL));
    cli_progress_set_name(pb, "Estimating ECDFs");
    nunprotect++;
  }

  for (int i=0; i < nr; i++) {
    SEXP          xR, uniqvR;
    int           nuniqv;
    double*       x;
    double*       uniqv;
    double*       ecdfuniqv;
    int           sum;
    double*       e1_p;
    const double* e2_p;
    int*          mt;
    int*          tab;
    int           nnas;

    if (verbose) { /* show progress */
      if (i % 100 == 0 && CLI_SHOULD_TICK)
        cli_progress_set(pb, i);
    }

    /* remove consecutive repeated elements */
    PROTECT(uniqvR = allocVector(REALSXP, nc));
    PROTECT(xR = allocVector(REALSXP, nc));
    uniqv = REAL(uniqvR);
    x = REAL(xR);
    nnas = 0; /* non-NA values */
    for (int j=0; j < nc; j++) {
#ifdef LONG_VECTOR_SUPPORT
      R_xlen_t idx = nr * j + i;
#else
      int idx = nr * j + i;
#endif
      x[j] = X[idx];
      if (!ISNA(x[j])) {
        uniqv[nnas] = x[j];
        nnas++;
      }
    }

    /* qsort(uniqv, nnas, sizeof(double), dbl_cmp); */
    R_qsort(uniqv, (size_t) 1, (size_t) nnas);
    e1_p = uniqv;
    e2_p = e1_p + 1;
    nuniqv = 0;
    if (nnas > 0)
      nuniqv = 1;
    for (int j=0; j < nnas-1; j++) { /* -1 for e2_p = e1_p + 1 */
      if (*e2_p != *e1_p) {
        *(++e1_p) = *e2_p;
        nuniqv++;
      }
      e2_p++;
    }

    /* match original values to sorted unique values */
    /* consider adding LONG_VECTOR_SUPPORT */
    mt = INTEGER(match_int(xR, uniqvR)); /* 1-based! */

    /* tabulate matches */
    /* consider adding LONG_VECTOR_SUPPORT */
    tab = R_Calloc(nuniqv, int); /* assuming zeroes are set */
    for (int j=0; j < nc; j++)
      if (mt[j] != NA_INTEGER && mt[j] > 0 && mt[j] <= nuniqv)
        tab[mt[j] - 1]++;

    /* cumulative sum to calculate ecdf values */
    /* consider adding LONG_VECTOR_SUPPORT */
    ecdfuniqv = R_Calloc(nuniqv, double); /* assuming zeroes are set */
    sum = 0;
    for (int j=0; j < nuniqv; j++) {
      sum = sum + tab[j];
      ecdfuniqv[j] = ((double) sum) / ((double) nc);
    }

    /* set ecdf values on the corresponding positions
     * of the output dense matrix */
    ecdf_vals = REAL(ecdfRobj);

    for (int j=0; j < nc; j++) {
#ifdef LONG_VECTOR_SUPPORT
      R_xlen_t idx = nr * j + i;
#else
      int idx = nr * j + i;
#endif
      if (!ISNA(X[idx]))
        ecdf_vals[idx] = ecdfuniqv[mt[j]-1];
      else
        ecdf_vals[idx] = NA_REAL;
    }

    R_Free(ecdfuniqv);
    R_Free(tab);

    UNPROTECT(2); /* xR uniqvR */
  }

  if (verbose)
    cli_progress_done(pb);

  UNPROTECT(nunprotect); /* XR ecdfRobj pb */

  return(ecdfRobj);
}
