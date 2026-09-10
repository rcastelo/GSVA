#include "ranks.h"

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
            SVT_SparseArray_svtSym,
            GSVA_attrNAsSym;


/* fetch integer column from a dense matrix XR of type integer
 * nr - number of rows
 * j - 0-based column to fetch
 * col - array where to store the column, assuming is initialized to zeroes
 * returned value - number of rows
 */
int
fetch_intcol_intmatrix(SEXP XR, int nr, int j, int* col) {
  int* X=INTEGER(XR);

  Memcpy(col, X+nr*j, (size_t) nr);

  return nr;
}

/* fetch integer column from a dense matrix XR of type double
 * nr - number of rows
 * j - 0-based column to fetch
 * col - array where to store the column, assuming is initialized to zeroes
 * returned value - number of rows
 */
int
fetch_intcol_dblmatrix(SEXP XR, int nr, int j, int* col) {
  double* X=REAL(XR);

  for (int i=0; i < nr; i++)
    col[i] = ISNAN(X[nr*j+i]) ? NA_INTEGER : (int) X[nr*j+i];

  return nr;
}

/* fetch integer column from a sparse dgCMatrix XR (always of type double)
 * nr - number of rows
 * j - 0-based column to fetch
 * col - array where to store the column, assuming is initialized to zeroes
 * returned value - number of nonzero values
 */
int
fetch_intcol_dgCMatrix(SEXP XCspR, int nr, int j, int* col) {
  int*    XCsp_i;
  int*    XCsp_p;
  double* XCsp_x;

  XCsp_i = INTEGER(GET_SLOT(XCspR, Matrix_iSym));
  XCsp_p = INTEGER(GET_SLOT(XCspR, Matrix_pSym));
  XCsp_x = REAL(GET_SLOT(XCspR, Matrix_xSym));

  /* put the sparse column into a dense vector */
  for (int i=XCsp_p[j]; i < XCsp_p[j+1]; i++)
    col[XCsp_i[i]] = ISNAN(XCsp_x[i]) ? NA_INTEGER : (int) XCsp_x[i];

  return XCsp_p[j+1]-XCsp_p[j];
}

/* fetch integer column from a sparse SVT_SparseMatrix XR of type integer
 * nr - number of rows
 * j - 0-based column to fetch
 * col - array where to store the column, assuming is initialized to zeroes
 * returned value - number of nonzero values
 */
int
fetch_intcol_intSVT_SparseMatrix(SEXP XsvtR, int nr, int j, int* col) {
  SEXP Xsvt_SVT;
  SEXP svtLeaf;
  int  nnz = 0;

  Xsvt_SVT = GET_SLOT(XsvtR, SVT_SparseArray_svtSym);
  svtLeaf = VECTOR_ELT(Xsvt_SVT, j);

  /* put the sparse column into a dense vector */
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
        col[offsets[i]] = vals[i];
    } else { /* lacunar */
      for (int i=0; i < noffsets; i++)
        col[offsets[i]] = 1;
    }
    nnz = noffsets;
  }

  return nnz;
}

/* fetch integer column from a sparse SVT_SparseMatrix XR of type double
 * nr - number of rows
 * j - 0-based column to fetch
 * col - array where to store the column, assuming is initialized to zeroes
 * returned value - number of nonzero values
 */
int
fetch_intcol_dblSVT_SparseMatrix(SEXP XsvtR, int nr, int j, int* col) {
  SEXP Xsvt_SVT;
  SEXP svtLeaf;
  int  nnz = 0;

  Xsvt_SVT = GET_SLOT(XsvtR, SVT_SparseArray_svtSym);
  svtLeaf = VECTOR_ELT(Xsvt_SVT, j);

  /* put the sparse column into a dense vector */
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
        col[offsets[i]] = ISNAN(vals[i]) ? NA_INTEGER : (int) vals[i];
    } else { /* lacunar */
      for (int i=0; i < noffsets; i++)
        col[offsets[i]] = 1;
    }
    nnz = noffsets;
  }

  return nnz;
}

FetchColFunDef
find_dim_and_fetchcolfun(SEXP XR, Rboolean intrnks, int** dim) {
  FetchColFunDef fetch_col;
  SEXP        classR = eval(lang2(install("class"), XR), R_BaseEnv);
  const char* class = CHAR(STRING_ELT(classR, 0));

  if (!strcmp(class, "matrix")) {
    fetch_col = intrnks ? &fetch_intcol_intmatrix : &fetch_intcol_dblmatrix;
    *dim = INTEGER(getAttrib(XR, R_DimSymbol));
  } else if (!strcmp(class, "dgCMatrix")) {
    fetch_col = &fetch_intcol_dgCMatrix;
    *dim = INTEGER(GET_SLOT(XR, Matrix_DimSym));
  } else if (!strcmp(class, "SVT_SparseMatrix")) {
    fetch_col = intrnks ? &fetch_intcol_intSVT_SparseMatrix :
                          &fetch_intcol_dblSVT_SparseMatrix;
    *dim = INTEGER(GET_SLOT(XR, SVT_SparseArray_dimSym));
  } else
    error("input class %s cannot be handled yet.", class);

  return fetch_col;
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

/* j is a 0-based column index on ranksR
 * returned value - 1 if some rank is out of range, 2 if some rank is missing,
 *                  which is not expected here, and 0 otherwise
 */
int
ranks2stats(SEXP ranksR, int p, int n, int j, Rboolean sparse,
            FetchColFunDef fetch_col,
            int* decordstat_col, double* symrnkstat_col,
            int* nzeros_col, double* zerosymrnkstat) {
  int* r = R_Calloc(p, int);       /* assume 0s are set */
  int* r_dense = R_Calloc(p, int); /* assume 0s are set */
  int  nnz, nzs;
  int  irankflag=0;

  nnz = (*fetch_col)(ranksR, p, j, r);
  nzs = p - nnz;

  if (nzs > 0) { /* if ranks have zeros, then input is a sparse matrix */
    int  k = 1;

    for (int i=0; i < p && !irankflag; i++) {
      /* missing ranks are not expected when calling this function */
      irankflag = r[i] == NA_INTEGER ? 2 : invalid_rank(r[i], p, nnz, 0);
      if (!irankflag) {
        if (r[i] == 0)       /* sparse ranks into dense ranks */
          r_dense[i] = k++;
        else
          r_dense[i] = r[i] + nzs;
      }
    }
  } else         /* input is a dense matrix */
    for (int i=0; i < p && !irankflag; i++) {
      /* missing ranks are not expected when calling this function */
      irankflag = r[i] == NA_INTEGER ? 2 : invalid_rank(r[i], p, nnz, 0);
      if (!irankflag)
        r_dense[i] = r[i];
    }

  if (irankflag) {
    R_Free(r_dense);
    R_Free(r);
    return irankflag;
  }

  /* dense ranks into decreasing order statistics */
  for (int i=0; i < p; i++)
    decordstat_col[i] = p - r_dense[i] + 1;

  *nzeros_col = *zerosymrnkstat = 0;
  if (nzs > 0 && sparse) {
    for (int i=0; i < p; i++) {
      double nnz1div2 = ((double) (nnz+1)) / 2.0;  /* nnz is in fact max(r)   */
      if (r[i] == 0) {                             /* in sparse regime zeros  */
        symrnkstat_col[i] = fabs(nnz1div2 - 1.0);  /* get same sym rank stat  */
        (*nzeros_col)++;                           /* count zeros  */
        *zerosymrnkstat = symrnkstat_col[i];       /* store zero sym rank stat*/
      } else                                       /* nonzero ranks shift one */
        symrnkstat_col[i] = fabs(nnz1div2 - (double) (r[i] + 1));
    }
  } else {
    for (int i=0; i < p; i++)
      symrnkstat_col[i] = fabs(((double) p) / 2.0 - ((double) r_dense[i]));
  }

  R_Free(r_dense);
  R_Free(r);

  return 0;
}

/* j is a 0-based column index on ranksR
 * regular reminder that ISNA() only applies to numeric values of
 * type double and missingess of integer values should be tested
 * by equality to NA_INTEGER, see
 * https://cran.r-project.org/doc/manuals/r-devel/R-exts.html#Missing-and-special-values-1
 * returned value - 1 if some rank is out of range, 0 otherwise
 */
int
ranks2stats_nas(SEXP ranksR, int p, int n, int j, Rboolean sparse,
                FetchColFunDef fetch_col,
                int* decordstat_col, double* symrnkstat_col,
                int* nzeros_col, double* zerosymrnkstat, int* nnas) {
  int* r = R_Calloc(p, int);       /* assume 0s are set */
  int* r_dense = R_Calloc(p, int); /* assume 0s are set */
  int  nnz, nzs;
  Rboolean irankflag=FALSE;

  nnz = (*fetch_col)(ranksR, p, j, r);
  *nnas = 0;
  for (int i=0; i < p; i++)
    if (r[i] == NA_INTEGER)
      (*nnas)++;

  nzs = p - nnz;

  if (nzs > 0) { /* if ranks have zeroes, then input is a sparse matrix */
    int  k = 1;

    for (int i=0; i < p && !irankflag; i++) {
      irankflag = invalid_rank(r[i], p, nnz, *nnas);
      if (!irankflag) {
        if (r[i] != NA_INTEGER) {
          if (r[i] == 0)       /* sparse ranks into dense ranks */
            r_dense[i] = k++;
          else
            r_dense[i] = r[i] + nzs;
        } else
          r_dense[i] = NA_INTEGER;
      }
    }
  } else         /* input is a dense matrix */
    for (int i=0; i < p && !irankflag; i++) {
      irankflag = invalid_rank(r[i], p, nnz, *nnas);
      if (!irankflag)
        r_dense[i] = r[i];
    }

  if (irankflag) {
    R_Free(r_dense);
    R_Free(r);
    return 1;
  }

  /* dense ranks into decreasing order statistics */
  for (int i=0; i < p; i++)
    decordstat_col[i] = r_dense[i] == NA_INTEGER ? NA_INTEGER : p - *nnas - r_dense[i] + 1;

  *nzeros_col = *zerosymrnkstat = 0;
  if (nzs > 0 && sparse) {
    for (int i=0; i < p; i++) {
      double nnz1div2 = ((double) (nnz-*nnas+1)) / 2.0; /* nnz is in fact max(r) */
      if (r[i] != NA_INTEGER) {
        if (r[i] == 0) {                             /* in sparse regime zeroes */
          symrnkstat_col[i] = fabs(nnz1div2 - 1.0);  /* get same sym rank stat  */
          (*nzeros_col)++;                           /* count zeros             */
          *zerosymrnkstat = symrnkstat_col[i];       /* store zero sym rank stat*/
        } else                                       /* nonzero ranks shift one */
          symrnkstat_col[i] = fabs(nnz1div2 - (double) (r[i] + 1));
      } else
          symrnkstat_col[i] = NA_REAL;
    }
  } else {
    for (int i=0; i < p; i++)
      if (r_dense[i] != NA_INTEGER)
        symrnkstat_col[i] = fabs(((double) (p - *nnas)) / 2.0 - ((double) r_dense[i]));
      else
        symrnkstat_col[i] = NA_REAL;
  }

  R_Free(r_dense);
  R_Free(r);

  return 0;
}
