#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <cli/progress.h>

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

int
fetch_row_nzvals(SEXP svtR, int i, int itypevals, int* whimin1,
                 int* inzvals, double* dnzvals, int* nzcols);

SEXP
row_rngs_nzrngs_RsparseMatrix_R(SEXP XRspR, SEXP verboseR) {
  SEXP     rngR;
  int*     XRsp_dim;
  int*     XRsp_p;
  Rboolean verbose=asLogical(verboseR);
  double*  XRsp_x;
  int      nr, nc;
  SEXP     pb=R_NilValue;
  double*  rng;
  int      nunprotect=0;

  PROTECT(XRspR); nunprotect++;

  XRsp_dim = INTEGER(GET_SLOT(XRspR, Matrix_DimSym));
  nr = XRsp_dim[0]; /* number of rows */
  nc = XRsp_dim[1]; /* number of columns */

  XRsp_p = INTEGER(GET_SLOT(XRspR, Matrix_pSym));
  XRsp_x = REAL(GET_SLOT(XRspR, Matrix_xSym));

  PROTECT(rngR = allocMatrix(REALSXP, nr, 4)); nunprotect++;

  if (verbose) {
    pb = PROTECT(cli_progress_bar(nr, NULL));
    cli_progress_set_name(pb, "Exploring rows");
    nunprotect++;
  }

  rng = REAL(rngR);

  for (int i=0; i < nr; i++) {
    int    nv;
    double min=NA_REAL;
    double max=NA_REAL;
    double nzmin=NA_REAL;
    double nzmax=NA_REAL;

    if (verbose) { /* show progress */
      if (i % 100 == 0 && CLI_SHOULD_TICK)
        cli_progress_set(pb, i);
    }

    /* number of nonzero values in the i-th row */
    nv = XRsp_p[i+1]-XRsp_p[i];
    if (nv < nc)     /* if there is at least one zero in the row */
      min = max = 0; /* then we can use it to initialize min and max */

    for (int j=XRsp_p[i]; j < XRsp_p[i+1]; j++) {
      if (!ISNA(XRsp_x[j])) {
        if (ISNA(max) || XRsp_x[j] > max)
          max = XRsp_x[j];
        if (ISNA(min) || XRsp_x[j] < min)
          min = XRsp_x[j];
        if (ISNA(nzmax) || XRsp_x[j] > nzmax)
          nzmax = XRsp_x[j];
        if (ISNA(nzmin) || XRsp_x[j] < nzmin)
          nzmin = XRsp_x[j];
      }
    }

    rng[i] = min;
    rng[nr + i] = max;
    rng[nr * 2 + i] = nzmin;
    rng[nr * 3 + i] = nzmax;
  }

  if (verbose)
    cli_progress_done(pb);

  UNPROTECT(nunprotect); /* XRspR rngR pb */

  return(rngR);
}

/* this function was only here for testing purposes */
/*
SEXP
row_rngs_nzrngs_SVT_SparseMatrix_R(SEXP XsvtR, SEXP verboseR) {
  SEXP        rngR;
  SEXP        Xsvt_dimR;
  const char* Xsvt_type;
  int*        Xsvt_dim;
  SEXP        Xsvt_SVT;
  int         itypevals;
  Rboolean    verbose=asLogical(verboseR);
  int         nr, nc;
  int*        whimin1;
  int*        inzvals=NULL;
  double*     dnzvals=NULL;
  SEXP        pb=R_NilValue;
  double*     rng;
  int         nunprotect=0;

  PROTECT(XsvtR); nunprotect++;

  Xsvt_dimR = GET_SLOT(XsvtR, SVT_SparseArray_dimSym);
  if (length(Xsvt_dimR) > 2)
    error("the input SVT_SparseMatrix object can only have two dimensions and has %d",
          length(Xsvt_dimR));

  Xsvt_dim = INTEGER(Xsvt_dimR);
  nr = Xsvt_dim[0];  * number of rows * 
  nc = Xsvt_dim[1];  * number of columns * 

  Xsvt_type = CHAR(STRING_ELT(getAttrib(XsvtR, SVT_SparseArray_typeSym), 0));
  Xsvt_SVT = GET_SLOT(XsvtR, SVT_SparseArray_svtSym);

  whimin1 = R_Calloc(nc, int);  * assuming values are set to 0s * 

  itypevals = 0;
  if (!strcmp(Xsvt_type, "integer")) {
    itypevals = 1;
    inzvals = R_Calloc(nc, int);
  } else
    dnzvals = R_Calloc(nc, double);

  PROTECT(rngR = allocMatrix(REALSXP, nr, 4)); nunprotect++;

  if (verbose) {
    pb = PROTECT(cli_progress_bar(nr, NULL));
    cli_progress_set_name(pb, "Exploring rows");
    nunprotect++;
  }

  rng = REAL(rngR);

  for (int i=0; i < nr; i++) {
    int    nv;
    double min=NA_REAL;
    double max=NA_REAL;
    double nzmin=NA_REAL;
    double nzmax=NA_REAL;

    if (verbose) {  * show progress * 
      if (i % 100 == 0 && CLI_SHOULD_TICK)
        cli_progress_set(pb, i);
    }

     * fetch nonzero values in the i-th row * 
    nv = fetch_row_nzvals(Xsvt_SVT, i, itypevals, whimin1,
                          inzvals, dnzvals, NULL);

    if (nv < nc)      * if there is at least one zero in the row * 
      min = max = 0;  * then we can use it to initialize min and max * 
     
    for (int j=0; j < nv; j++) {
      double x = itypevals ? ((double) inzvals[j]) : dnzvals[j];
      if (!ISNA(x)) {
        if (ISNA(max) || x > max)
          max = x;
        if (ISNA(min) || x < min)
          min = x;
        if (ISNA(nzmax) || x > nzmax)
          nzmax = x;
        if (ISNA(nzmin) || x < nzmin)
          nzmin = x;
      }
    }

    rng[i] = min;
    rng[nr + i] = max;
    rng[nr * 2 + i] = nzmin;
    rng[nr * 3 + i] = nzmax;
  }

  R_Free(whimin1);
  if (itypevals)
    R_Free(inzvals);
  else
    R_Free(dnzvals);

  if (verbose)
    cli_progress_done(pb);

  UNPROTECT(nunprotect);  * XsvtR rngR pb * 

  return(rngR);
}
*/

/* this function was only here for testing purposes */
/*
SEXP
col_rngs_nzrngs_SVT_SparseMatrix_R(SEXP XsvtR, SEXP verboseR) {
  SEXP        rngR;
  SEXP        Xsvt_dimR;
  const char* Xsvt_type;
  int*        Xsvt_dim;
  SEXP        Xsvt_SVT;
  int         itypevals;
  Rboolean    verbose=asLogical(verboseR);
  int         nr, nc;
  SEXP        pb=R_NilValue;
  double*     rng;
  int         nunprotect=0;

  PROTECT(XsvtR); nunprotect++;

  Xsvt_dimR = GET_SLOT(XsvtR, SVT_SparseArray_dimSym);
  if (length(Xsvt_dimR) > 2)
    error("the input SVT_SparseMatrix object can only have two dimensions and has %d",
          length(Xsvt_dimR));

  Xsvt_dim = INTEGER(Xsvt_dimR);
  nr = Xsvt_dim[0];  * number of rows * 
  nc = Xsvt_dim[1];  * number of columns * 

  Xsvt_type = CHAR(STRING_ELT(getAttrib(XsvtR, SVT_SparseArray_typeSym), 0));
  Xsvt_SVT = GET_SLOT(XsvtR, SVT_SparseArray_svtSym);

  itypevals = 0;
  if (!strcmp(Xsvt_type, "integer")) {
    itypevals = 1;
  }

  PROTECT(rngR = allocMatrix(REALSXP, nc, 4)); nunprotect++;

  if (verbose) {
    pb = PROTECT(cli_progress_bar(nr, NULL));
    cli_progress_set_name(pb, "Exploring rows");
    nunprotect++;
  }

  rng = REAL(rngR);

  for (int i=0; i < nc; i++) {
    SEXP   svtLeaf;
    double min=NA_REAL;
    double max=NA_REAL;
    double nzmin=NA_REAL;
    double nzmax=NA_REAL;

    if (verbose) {  * show progress * 
      if (i % 100 == 0 && CLI_SHOULD_TICK)
        cli_progress_set(pb, i);
    }

    if (length(Xsvt_SVT) > 0) {  * matrix not empty * 
      svtLeaf = VECTOR_ELT(Xsvt_SVT, i);
      if (svtLeaf != R_NilValue) {
        SEXP    valsR = VECTOR_ELT(svtLeaf, 0);
        SEXP    offsetsR = VECTOR_ELT(svtLeaf, 1);
        int     nvals = length(valsR);
        int     noffsets = length(offsetsR);
        int*    inzvals;
        double* dnzvals;

        if (noffsets < nr)      * if there is at least one zero in the column * 
          min = max = 0;        * then we can use it to initialize min and max * 

        if (nvals > 0) {
          if (itypevals)
            inzvals = INTEGER(valsR);
          else
            dnzvals = REAL(valsR);

          for (int j=0; j < nvals; j++) {
            double x = itypevals ? ((double) inzvals[j]) : dnzvals[j];
            if (!ISNA(x)) {
              if (ISNA(max) || x > max)
                max = x;
              if (ISNA(min) || x < min)
                min = x;
              if (ISNA(nzmax) || x > nzmax)
                nzmax = x;
              if (ISNA(nzmin) || x < nzmin)
                nzmin = x;
            }
          }
        } else {
          max = nzmin = nzmax = 1;
          if (ISNA(min))
            min = 1;
        }
      } else
        min = max = 0;
    } else
      min = max = 0;

    rng[i] = min;
    rng[nc + i] = max;
    rng[nc * 2 + i] = nzmin;
    rng[nc * 3 + i] = nzmax;
  }

  if (verbose)
    cli_progress_done(pb);

  UNPROTECT(nunprotect);  * XsvtR rngR pb * 

  return(rngR);
}
*/

SEXP
rowbycols_rngs_nzrngs_SVT_SparseMatrix_R(SEXP XsvtR, SEXP verboseR) {
  SEXP        rngR;
  SEXP        Xsvt_dimR;
  const char* Xsvt_type;
  int*        Xsvt_dim;
  SEXP        Xsvt_SVT;
  int         itypevals;
  Rboolean    verbose=asLogical(verboseR);
  int         nr, nc;
  SEXP        pb=R_NilValue;
  double*     rng;
  int         nunprotect=0;

  PROTECT(XsvtR); nunprotect++;

  Xsvt_dimR = GET_SLOT(XsvtR, SVT_SparseArray_dimSym);
  if (length(Xsvt_dimR) > 2)
    error("the input SVT_SparseMatrix object can only have two dimensions and has %d",
          length(Xsvt_dimR));

  Xsvt_dim = INTEGER(Xsvt_dimR);
  nr = Xsvt_dim[0]; /* number of rows */
  nc = Xsvt_dim[1]; /* number of columns */

  Xsvt_type = CHAR(STRING_ELT(getAttrib(XsvtR, SVT_SparseArray_typeSym), 0));
  Xsvt_SVT = GET_SLOT(XsvtR, SVT_SparseArray_svtSym);

  itypevals = 0;
  if (!strcmp(Xsvt_type, "integer")) {
    itypevals = 1;
  }

  PROTECT(rngR = allocMatrix(REALSXP, nr, 4)); nunprotect++;

  if (verbose) {
    pb = PROTECT(cli_progress_bar(nr, NULL)); nunprotect++;
    cli_progress_set_name(pb, "Exploring rows");
  }

  rng = REAL(rngR);

  for (int i=0; i < nr; i++) {
    rng[i] = rng[nr + i] = 0;                    /* init min max */
    rng[nr * 2 + i] = rng[nr * 3 + i] = NA_REAL; /* init nzmin nzmax */
  }

  if (length(Xsvt_SVT) == 0) { /* input matrix is empty */
    if (verbose)
      cli_progress_done(pb);

    UNPROTECT(nunprotect); /* XsvtR rngR pb */

    return(rngR);
  }

  for (int i=0; i < nc; i++) {
    SEXP svtLeaf = VECTOR_ELT(Xsvt_SVT, i);

    if (verbose) { /* show progress */
      if (i % 100 == 0 && CLI_SHOULD_TICK)
        cli_progress_set(pb, i);
    }

    if (svtLeaf != R_NilValue) {
      SEXP    valsR = VECTOR_ELT(svtLeaf, 0);
      SEXP    offsetsR = VECTOR_ELT(svtLeaf, 1);
      int*    offsets = INTEGER(offsetsR);
      int     nvals = length(valsR);
      int     noffsets = length(offsetsR);
      int*    inzvals;
      double* dnzvals;

      if (nvals > 0) {
        if (itypevals)
          inzvals = INTEGER(valsR);
        else
          dnzvals = REAL(valsR);

        for (int j=0; j < nvals; j++) {
          double x = itypevals ? ((double) inzvals[j]) : dnzvals[j];
          double min = rng[offsets[j]];
          double max = rng[nr + offsets[j]];
          double nzmin = rng[nr * 2 + offsets[j]];
          double nzmax = rng[nr * 3 + offsets[j]];

          if (!ISNA(x)) {
            if (ISNA(min) || x < min)
              rng[offsets[j]] = x;
            if (ISNA(max) || x > max)
              rng[nr + offsets[j]] = x;
            if (ISNA(nzmin) || x < nzmin)
              rng[nr * 2 + offsets[j]] = x;
            if (ISNA(nzmax) || x > nzmax)
              rng[nr * 3 + offsets[j]] = x;
          }
        }
      } else { /* lacunar leaf */
        for (int j=0; j < noffsets; j++) {
          double min = rng[offsets[j]];
          double max = rng[nr + offsets[j]];
          double nzmin = rng[nr * 2 + offsets[j]];
          double nzmax = rng[nr * 3 + offsets[j]];
          if (ISNA(min) || min > 1)
            rng[offsets[j]] = 1;
          if (ISNA(max) || max < 1)
            rng[nr + offsets[j]] = 1;
          if (ISNA(nzmin) || min > 1)
            rng[nr * 2 + offsets[j]] = 1;
          if (ISNA(nzmax) || max < 1)
            rng[nr * 3 + offsets[j]] = 1;
        }
      }
    }
  }

  if (verbose)
    cli_progress_done(pb);

  UNPROTECT(nunprotect); /* XsvtR rngR pb */

  return(rngR);
}
