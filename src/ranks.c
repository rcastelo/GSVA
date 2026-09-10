#include "ranks.h"
#include <string.h>

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


/* fetch integer column from a dense matrix of type integer */
static int
fetch_intcol_intmatrix(struct ranks_ctx_s* ctx, int j, int* col) {
  int* X = (int*)ctx->u.data;
  Memcpy(col, X + ctx->p * j, (size_t)ctx->p);
  return ctx->p;
}

/* fetch integer column from a dense matrix of type double */
static int
fetch_intcol_dblmatrix(struct ranks_ctx_s* ctx, int j, int* col) {
  double* X = (double*)ctx->u.data;
  for (int i = 0; i < ctx->p; i++)
    col[i] = ISNAN(X[ctx->p * j + i]) ? NA_INTEGER : (int)X[ctx->p * j + i];
  return ctx->p;
}

/* fetch integer column from a sparse dgCMatrix */
static int
fetch_intcol_dgCMatrix(struct ranks_ctx_s* ctx, int j, int* col) {
  const int*    i  = ctx->u.dgc.i;
  const int*    p  = ctx->u.dgc.p;
  const double* x  = ctx->u.dgc.x;

  memset(col, 0, (size_t)ctx->p * sizeof(int));
  for (int idx = p[j]; idx < p[j + 1]; idx++)
    col[i[idx]] = ISNAN(x[idx]) ? NA_INTEGER : (int)x[idx];

  return p[j + 1] - p[j];
}

/* fetch integer column from a sparse SVT_SparseMatrix of type integer */
static int
fetch_intcol_intSVT_SparseMatrix(struct ranks_ctx_s* ctx, int j, int* col) {
  SEXP Xsvt_SVT = ctx->u.svt;
  SEXP svtLeaf  = VECTOR_ELT(Xsvt_SVT, j);
  int  nnz      = 0;

  memset(col, 0, (size_t)ctx->p * sizeof(int));
  if (svtLeaf != R_NilValue) {
    SEXP valsR        = VECTOR_ELT(svtLeaf, 0);
    SEXP offsetsR     = VECTOR_ELT(svtLeaf, 1);
    int  nvals        = length(valsR);
    int  noffsets     = length(offsetsR);
    int* vals;
    int* offsets = INTEGER(offsetsR);

    if (nvals > 0) {
      vals = INTEGER(valsR);
      for (int i = 0; i < nvals; i++)
        col[offsets[i]] = vals[i];
    } else { /* lacunar */
      for (int i = 0; i < noffsets; i++)
        col[offsets[i]] = 1;
    }
    nnz = noffsets;
  }

  return nnz;
}

/* fetch integer column from a sparse SVT_SparseMatrix of type double */
static int
fetch_intcol_dblSVT_SparseMatrix(struct ranks_ctx_s* ctx, int j, int* col) {
  SEXP  Xsvt_SVT = ctx->u.svt;
  SEXP  svtLeaf  = VECTOR_ELT(Xsvt_SVT, j);
  int   nnz      = 0;

  memset(col, 0, (size_t)ctx->p * sizeof(int));
  if (svtLeaf != R_NilValue) {
    SEXP   valsR        = VECTOR_ELT(svtLeaf, 0);
    SEXP   offsetsR     = VECTOR_ELT(svtLeaf, 1);
    int    nvals        = length(valsR);
    int    noffsets     = length(offsetsR);
    double* vals;
    int*    offsets = INTEGER(offsetsR);

    if (nvals > 0) {
      vals = REAL(valsR);
      for (int i = 0; i < nvals; i++)
        col[offsets[i]] = ISNAN(vals[i]) ? NA_INTEGER : (int)vals[i];
    } else { /* lacunar */
      for (int i = 0; i < noffsets; i++)
        col[offsets[i]] = 1;
    }
    nnz = noffsets;
  }

  return nnz;
}

ranks_ctx_t*
ranks_ctx_create(SEXP XR, Rboolean intrnks, Rboolean sparse) {
  ranks_ctx_t* ctx     = (ranks_ctx_t*)R_alloc(1, sizeof(ranks_ctx_t));
  SEXP         classR  = eval(lang2(install("class"), XR), R_BaseEnv);
  const char*  class   = CHAR(STRING_ELT(classR, 0));
  int*         dim;

  ctx->sparse = sparse;

  if (!strcmp(class, "matrix")) {
    dim = INTEGER(getAttrib(XR, R_DimSymbol));
    if (intrnks) {
      ctx->type      = RANKSTYPE_MATRIX_INT;
      ctx->u.data    = (const void*)INTEGER(XR);
      ctx->fetch_col = &fetch_intcol_intmatrix;
    } else {
      ctx->type      = RANKSTYPE_MATRIX_DBL;
      ctx->u.data    = (const void*)REAL(XR);
      ctx->fetch_col = &fetch_intcol_dblmatrix;
    }
  } else if (!strcmp(class, "dgCMatrix")) {
    dim               = INTEGER(GET_SLOT(XR, Matrix_DimSym));
    ctx->type         = RANKSTYPE_DGC;
    ctx->u.dgc.i      = INTEGER(GET_SLOT(XR, Matrix_iSym));
    ctx->u.dgc.p      = INTEGER(GET_SLOT(XR, Matrix_pSym));
    ctx->u.dgc.x      = REAL(GET_SLOT(XR, Matrix_xSym));
    ctx->fetch_col    = &fetch_intcol_dgCMatrix;
  } else if (!strcmp(class, "SVT_SparseMatrix")) {
    dim                = INTEGER(GET_SLOT(XR, SVT_SparseArray_dimSym));
    ctx->u.svt         = GET_SLOT(XR, SVT_SparseArray_svtSym);
    if (intrnks) {
      ctx->type        = RANKSTYPE_SVT_INT;
      ctx->fetch_col   = &fetch_intcol_intSVT_SparseMatrix;
    } else {
      ctx->type        = RANKSTYPE_SVT_DBL;
      ctx->fetch_col   = &fetch_intcol_dblSVT_SparseMatrix;
    }
  } else {
    error("input class %s cannot be handled yet.", class);
  }

  ctx->p    = dim[0];
  ctx->n    = dim[1];

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

int
ranks2stats(ranks_ctx_t* ctx, int j, int* decordstat, double* symrnkstat,
            int* nzeros_col, double* zerosymrnkstat) {
  int* r         = (int*)R_alloc(ctx->p, sizeof(int));
  int* r_dense   = (int*)R_alloc(ctx->p, sizeof(int));
  int  nnz, nzs;
  int irankflag = 0;

  nnz = ctx->fetch_col(ctx, j, r);
  nzs = ctx->p - nnz;

  if (nzs > 0) { /* input is a sparse matrix */
    int k = 1;
    for (int i = 0; i < ctx->p && !irankflag; i++) {
      /* missing ranks are not expected when calling this function */
      irankflag = r[i] == NA_INTEGER ? 2 : invalid_rank(r[i], ctx->p, nnz, 0);
      if (!irankflag) {
        if (r[i] == 0)
          r_dense[i] = k++;
        else
          r_dense[i] = r[i] + nzs;
      }
    }
  } else { /* input is a dense matrix */
    for (int i = 0; i < ctx->p && !irankflag; i++) {
      /* missing ranks are not expected when calling this function */
      irankflag = r[i] == NA_INTEGER ? 2 : invalid_rank(r[i], ctx->p, nnz, 0);
      if (!irankflag)
        r_dense[i] = r[i];
    }
  }

  if (irankflag)
    return irankflag;

  /* dense ranks into decreasing order statistics */
  for (int i = 0; i < ctx->p; i++)
    decordstat[i] = ctx->p - r_dense[i] + 1;

  *nzeros_col      = 0;
  *zerosymrnkstat  = 0.0;
  if (nzs > 0 && ctx->sparse) {
    double nnz1div2 = ((double)(nnz + 1)) / 2.0;
    for (int i = 0; i < ctx->p; i++) {
      if (r[i] == 0) {
        symrnkstat[i] = fabs(nnz1div2 - 1.0);
        (*nzeros_col)++;
        *zerosymrnkstat = symrnkstat[i];
      }
      else
        symrnkstat[i] = fabs(nnz1div2 - (double)(r[i] + 1));
    }
  } else {
    for (int i = 0; i < ctx->p; i++)
      symrnkstat[i] = fabs(((double)ctx->p) / 2.0 - ((double)r_dense[i]));
  }
  return 0;
}

int
ranks2stats_nas(ranks_ctx_t* ctx, int j, int* decordstat, double* symrnkstat,
                int* nzeros_col, double* zerosymrnkstat, int* nnas) {
  int* r         = (int*)R_alloc(ctx->p, sizeof(int));
  int* r_dense   = (int*)R_alloc(ctx->p, sizeof(int));
  int  nnz, nzs;
  Rboolean irankflag = FALSE;

  nnz = ctx->fetch_col(ctx, j, r);
  *nnas = 0;
  for (int i = 0; i < ctx->p; i++)
    if (r[i] == NA_INTEGER)
      (*nnas)++;

  nzs = ctx->p - nnz;

  if (nzs > 0) { /* input is a sparse matrix */
    int k = 1;
    for (int i = 0; i < ctx->p && !irankflag; i++) {
      irankflag = invalid_rank(r[i], ctx->p, nnz, *nnas);
      if (!irankflag) {
        if (r[i] != NA_INTEGER) {
          if (r[i] == 0)
            r_dense[i] = k++;
          else
            r_dense[i] = r[i] + nzs;
        } else
          r_dense[i] = NA_INTEGER;
      }
    }
  } else { /* input is a dense matrix */
    for (int i = 0; i < ctx->p && !irankflag; i++) {
      irankflag = invalid_rank(r[i], ctx->p, nnz, *nnas);
      if (!irankflag)
        r_dense[i] = r[i];
    }
  }

  if (irankflag)
    return 1;

  /* dense ranks into decreasing order statistics */
  for (int i = 0; i < ctx->p; i++)
    decordstat[i] = r_dense[i] == NA_INTEGER ? NA_INTEGER
                                             : ctx->p - *nnas - r_dense[i] + 1;

  *nzeros_col      = 0;
  *zerosymrnkstat  = 0.0;
  if (nzs > 0 && ctx->sparse) {
    double nnz1div2 = ((double)(nnz - *nnas + 1)) / 2.0;
    for (int i = 0; i < ctx->p; i++) {
      if (r[i] != NA_INTEGER) {
        if (r[i] == 0) {
          symrnkstat[i] = fabs(nnz1div2 - 1.0);
          (*nzeros_col)++;
          *zerosymrnkstat = symrnkstat[i];
        }
        else
          symrnkstat[i] = fabs(nnz1div2 - (double)(r[i] + 1));
      }
      else
        symrnkstat[i] = NA_REAL;
    }
  } else {
    for (int i = 0; i < ctx->p; i++) {
      if (r_dense[i] != NA_INTEGER)
        symrnkstat[i] = fabs(((double)(ctx->p - *nnas)) / 2.0
                                  - ((double)r_dense[i]));
      else
        symrnkstat[i] = NA_REAL;
    }
  }
  return 0;
}
