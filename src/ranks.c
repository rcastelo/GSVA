#include <string.h>
#include <math.h>
#include "ranks.h"

/* global variables */
extern SEXP Matrix_DimSym,
            Matrix_xSym,
            Matrix_iSym,
            Matrix_pSym,
            SVT_SparseArray_dimSym,
            SVT_SparseArray_svtSym;

/* fetch integer column from a dense matrix of type integer
 * j - 0-based column to fetch
 * returned value - number of rows
 */
static int
fetch_intcol_intmatrix(ranks_ctx_t* ctx, int j) {
  const int* X = (const int*) ctx->u.data;
#ifdef LONG_VECTOR_SUPPORT
  R_xlen_t   offset = (R_xlen_t) ctx->p * j;
#else
  int        offset = ctx->p * j;
#endif

  Memcpy(ctx->r, X+offset, (size_t) ctx->p);

  return ctx->p;
}

/* fetch integer column from a dense matrix of type double
 * j - 0-based column to fetch
 * returned value - number of rows
 */
static int
fetch_intcol_dblmatrix(ranks_ctx_t* ctx, int j) {
  const double* X = (const double*) ctx->u.data;
#ifdef LONG_VECTOR_SUPPORT
  R_xlen_t      offset = (R_xlen_t) ctx->p * j;
#else
  int           offset = ctx->p * j;
#endif

  for (int i=0; i < ctx->p; i++)
    ctx->r[i] = ISNAN(X[offset+i]) ? NA_INTEGER : (int) X[offset+i];

  return ctx->p;
}

/* fetch integer column from a sparse dgCMatrix (always of type double)
 * j - 0-based column to fetch
 * returned value - number of nonzero values
 */
static int
fetch_intcol_dgCMatrix(ranks_ctx_t* ctx, int j) {
  const int*    XCsp_i = ctx->u.dgc.i;
  const int*    XCsp_p = ctx->u.dgc.p;
  const double* XCsp_x = ctx->u.dgc.x;

  /* put the sparse column into a dense vector */
  Memzero(ctx->r, (size_t) ctx->p);
  for (int i=XCsp_p[j]; i < XCsp_p[j+1]; i++)
    ctx->r[XCsp_i[i]] = ISNAN(XCsp_x[i]) ? NA_INTEGER : (int) XCsp_x[i];

  return XCsp_p[j+1]-XCsp_p[j];
}

/* fetch integer column from a sparse SVT_SparseMatrix of type integer
 * j - 0-based column to fetch
 * returned value - number of nonzero values
 */
static int
fetch_intcol_intSVT_SparseMatrix(ranks_ctx_t* ctx, int j) {
  SEXP svtLeaf = VECTOR_ELT(ctx->u.svt, j);
  int  nnz = 0;

  /* put the sparse column into a dense vector */
  Memzero(ctx->r, (size_t) ctx->p);
  if (svtLeaf != R_NilValue) {
    SEXP valsR = VECTOR_ELT(svtLeaf, 0);
    SEXP offsetsR = VECTOR_ELT(svtLeaf, 1);
    int  nvals = length(valsR);
    int  noffsets = length(offsetsR);
    int* vals;
    int* offsets = INTEGER(offsetsR);

    if (nvals > 0) {
      vals = INTEGER(valsR);
      for (int i=0; i < nvals; i++)
        ctx->r[offsets[i]] = vals[i];
    } else { /* lacunar */
      for (int i=0; i < noffsets; i++)
        ctx->r[offsets[i]] = 1;
    }
    nnz = noffsets;
  }

  return nnz;
}

/* fetch integer column from a sparse SVT_SparseMatrix of type double
 * j - 0-based column to fetch
 * returned value - number of nonzero values
 */
static int
fetch_intcol_dblSVT_SparseMatrix(ranks_ctx_t* ctx, int j) {
  SEXP svtLeaf = VECTOR_ELT(ctx->u.svt, j);
  int  nnz = 0;

  /* put the sparse column into a dense vector */
  Memzero(ctx->r, (size_t) ctx->p);
  if (svtLeaf != R_NilValue) {
    SEXP    valsR = VECTOR_ELT(svtLeaf, 0);
    SEXP    offsetsR = VECTOR_ELT(svtLeaf, 1);
    int     nvals = length(valsR);
    int     noffsets = length(offsetsR);
    double* vals;
    int*    offsets = INTEGER(offsetsR);

    if (nvals > 0) {
      vals = REAL(valsR);
      for (int i=0; i < nvals; i++)
        ctx->r[offsets[i]] = ISNAN(vals[i]) ? NA_INTEGER : (int) vals[i];
    } else { /* lacunar */
      for (int i=0; i < noffsets; i++)
        ctx->r[offsets[i]] = 1;
    }
    nnz = noffsets;
  }

  return nnz;
}

/* figures out the class of the rank matrix 'XR', sets the column fetching
 * function and the pointers to its data, and allocates the column buffer.
 * memory is allocated with R_alloc(), which is released by R when the .Call()
 * returns, also on error()
 */
ranks_ctx_t*
ranks_ctx_create(SEXP XR, Rboolean intrnks, Rboolean sparse) {
  ranks_ctx_t* ctx = (ranks_ctx_t*) R_alloc(1, sizeof(ranks_ctx_t));
  const char*  class = get_class_name(XR);
  int*         dim;

  ctx->sparse = sparse;

  if (!strcmp(class, "matrix")) {
    dim = INTEGER(getAttrib(XR, R_DimSymbol));
    if (intrnks) {
      ctx->type = RANKSTYPE_MATRIX_INT;
      ctx->u.data = (const void*) INTEGER(XR);
      ctx->fetch_col = &fetch_intcol_intmatrix;
    } else {
      ctx->type = RANKSTYPE_MATRIX_DBL;
      ctx->u.data = (const void*) REAL(XR);
      ctx->fetch_col = &fetch_intcol_dblmatrix;
    }
  } else if (!strcmp(class, "dgCMatrix")) {
    dim = INTEGER(GET_SLOT(XR, Matrix_DimSym));
    ctx->type = RANKSTYPE_DGC;
    ctx->u.dgc.i = INTEGER(GET_SLOT(XR, Matrix_iSym));
    ctx->u.dgc.p = INTEGER(GET_SLOT(XR, Matrix_pSym));
    ctx->u.dgc.x = REAL(GET_SLOT(XR, Matrix_xSym));
    ctx->fetch_col = &fetch_intcol_dgCMatrix;
  } else if (!strcmp(class, "SVT_SparseMatrix")) {
    dim = INTEGER(GET_SLOT(XR, SVT_SparseArray_dimSym));
    ctx->u.svt = GET_SLOT(XR, SVT_SparseArray_svtSym);
    if (intrnks) {
      ctx->type = RANKSTYPE_SVT_INT;
      ctx->fetch_col = &fetch_intcol_intSVT_SparseMatrix;
    } else {
      ctx->type = RANKSTYPE_SVT_DBL;
      ctx->fetch_col = &fetch_intcol_dblSVT_SparseMatrix;
    }
  } else
    error("input class %s cannot be handled yet.", class);

  ctx->p = dim[0];
  ctx->n = dim[1];
  ctx->r = (int*) R_alloc(ctx->p, sizeof(int));

  return ctx;
}

/* check that a rank 'r' of a column with 'p' rows, of which 'nnz' are
 * stored values (all 'p' in a dense column) and 'nnas' are missing, is
 * either zero, only when there are zeros in a sparse column, or between 1
 * and the number of stored nonmissing values, which is required to convert
 * it into a decreasing order statistic between 1 and 'p' minus 'nnas'.
 * returned value - 1 if the rank is out of range, 0 otherwise
 */
static Rboolean
invalid_rank(int r, int p, int nnz, int nnas) {
  int maxrank = nnz - nnas;

  if (r != NA_INTEGER) {
    if (r == 0) {
      if (nnz == p)   /* a dense column cannot have zero ranks */
        return TRUE;
    } else if (r < 1 || r > maxrank)
      return TRUE;
  }

  return FALSE;
}

/* j is a 0-based column index on the rank matrix; the ranks of the column are
 * validated before being converted in place, in the buffer of the context,
 * into dense ranks
 * returned value - 1 if some rank is out of range, 2 if some rank is missing,
 *                  which is not expected here, and 0 otherwise
 */
int
ranks2stats(ranks_ctx_t* ctx, int j, int* decordstat, double* symrnkstat,
            int* nzeros_col, double* zerosymrnkstat) {
  int* r = ctx->r;
  int  p = ctx->p;
  int  nnz, nzs;
  int  irankflag=0;

  nnz = ctx->fetch_col(ctx, j);
  nzs = p - nnz;

  for (int i=0; i < p && !irankflag; i++)
    /* missing ranks are not expected when calling this function */
    irankflag = r[i] == NA_INTEGER ? 2 : invalid_rank(r[i], p, nnz, 0);

  if (irankflag)
    return irankflag;

  *nzeros_col = *zerosymrnkstat = 0;
  if (nzs > 0) { /* if ranks have zeros, then input is a sparse matrix */
    int k = 1;

    if (ctx->sparse) {
      double nnz1div2 = ((double) (nnz+1)) / 2.0; /* nnz is in fact max(r)    */
      for (int i=0; i < p; i++) {
        if (r[i] == 0) {                          /* in sparse regime zeros   */
          symrnkstat[i] = fabs(nnz1div2 - 1.0);   /* get same sym rank stat   */
          (*nzeros_col)++;                        /* count zeros              */
          *zerosymrnkstat = symrnkstat[i];        /* store zero sym rank stat */
        } else                                    /* nonzero ranks shift one  */
          symrnkstat[i] = fabs(nnz1div2 - (double) (r[i] + 1));
      }
    }

    for (int i=0; i < p; i++) /* sparse ranks into dense ranks */
      r[i] = r[i] == 0 ? k++ : r[i] + nzs;
  }

  /* dense ranks into decreasing order statistics */
  for (int i=0; i < p; i++)
    decordstat[i] = p - r[i] + 1;

  if (!(nzs > 0 && ctx->sparse))
    for (int i=0; i < p; i++)
      symrnkstat[i] = fabs(((double) p) / 2.0 - ((double) r[i]));

  return 0;
}

/* j is a 0-based column index on the rank matrix
 * regular reminder that ISNA() only applies to numeric values of
 * type double and missingess of integer values should be tested
 * by equality to NA_INTEGER, see
 * https://cran.r-project.org/doc/manuals/r-devel/R-exts.html#Missing-and-special-values-1
 * returned value - 1 if some rank is out of range, 0 otherwise
 */
int
ranks2stats_nas(ranks_ctx_t* ctx, int j, int* decordstat, double* symrnkstat,
                int* nzeros_col, double* zerosymrnkstat, int* nnas) {
  int* r = ctx->r;
  int  p = ctx->p;
  int  nnz, nzs;
  Rboolean irankflag=FALSE;

  nnz = ctx->fetch_col(ctx, j);
  *nnas = 0;
  for (int i=0; i < p; i++)
    if (r[i] == NA_INTEGER)
      (*nnas)++;

  nzs = p - nnz;

  for (int i=0; i < p && !irankflag; i++)
    irankflag = invalid_rank(r[i], p, nnz, *nnas);

  if (irankflag)
    return 1;

  *nzeros_col = *zerosymrnkstat = 0;
  if (nzs > 0) { /* if ranks have zeroes, then input is a sparse matrix */
    int k = 1;

    if (ctx->sparse) {
      double nnz1div2 = ((double) (nnz-*nnas+1)) / 2.0; /* nnz is in fact max(r) */
      for (int i=0; i < p; i++) {
        if (r[i] != NA_INTEGER) {
          if (r[i] == 0) {                        /* in sparse regime zeroes  */
            symrnkstat[i] = fabs(nnz1div2 - 1.0); /* get same sym rank stat   */
            (*nzeros_col)++;                      /* count zeros              */
            *zerosymrnkstat = symrnkstat[i];      /* store zero sym rank stat */
          } else                                  /* nonzero ranks shift one  */
            symrnkstat[i] = fabs(nnz1div2 - (double) (r[i] + 1));
        } else
          symrnkstat[i] = NA_REAL;
      }
    }

    for (int i=0; i < p; i++) /* sparse ranks into dense ranks */
      if (r[i] != NA_INTEGER)
        r[i] = r[i] == 0 ? k++ : r[i] + nzs;
  }

  /* dense ranks into decreasing order statistics */
  for (int i=0; i < p; i++)
    decordstat[i] = r[i] == NA_INTEGER ? NA_INTEGER : p - *nnas - r[i] + 1;

  if (!(nzs > 0 && ctx->sparse)) {
    for (int i=0; i < p; i++) {
      if (r[i] != NA_INTEGER)
        symrnkstat[i] = fabs(((double) (p - *nnas)) / 2.0 - ((double) r[i]));
      else
        symrnkstat[i] = NA_REAL;
    }
  }

  return 0;
}
