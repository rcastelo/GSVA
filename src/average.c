#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>
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
            SVT_SparseArray_svtSym,
            GSVA_attrNAsSym;

const char* get_class_name(SEXP obj);

/* fetch a column from a dense matrix of type double
 * XR - input matrix
 * nr - number of rows in XR
 * j - 0-based column to fetch
 * col - pointer to where the returned column should be stored
 * idx - pointer to where the returned column indices should be stored
 *       (NULL for dense matrices)
 * returned value - number of rows
 */
int
fetchnzcol_dblmatrix(SEXP XR, int nr, int j, double* col, int* idx) {
    double* X=REAL(XR);

    Memcpy(col, X + nr * j, (size_t) nr);
    idx = NULL;

    return nr;
}

/* fetch column nonzero values from a dgCMatrix object of type double
 * XCspR - input matrix of type dgCMatrix
 * nr - number of rows in XCspR
 * j - 0-based column to fetch
 * col - pointer to where the returned column nonzero values should be stored
 * idx - pointer to where the returned column nonzero indices should be stored
 * returned value - number of rows
 */
int
fetchnzcol_dgCMatrix(SEXP XCspR, int nr, int j, double* col, int* idx) {
  int*    XCsp_i;
  int*    XCsp_p;
  double* XCsp_x;
  int     nnz;

  XCsp_i = INTEGER(GET_SLOT(XCspR, Matrix_iSym));
  XCsp_p = INTEGER(GET_SLOT(XCspR, Matrix_pSym));
  XCsp_x = REAL(GET_SLOT(XCspR, Matrix_xSym));

  nnz = XCsp_p[j+1] - XCsp_p[j];
  Memcpy(col, XCsp_x + XCsp_p[j], (size_t) nnz);
  Memcpy(idx, XCsp_i + XCsp_p[j], (size_t) nnz);

  return nnz;
}

/* fetch pointer to a column from a SVT_SparseMatrix object of type double
 * XsvtR - input matrix of type SVT_SparseMatrix
 * nr - number of rows in XsvtR
 * j - 0-based column to fetch
 * col - pointer to the returned column
 * idx - pointer to the returned column indices (NULL for dense matrix)
 * returned value - number of rows
 */
int
fetchnzcol_dblSVT_SparseMatrix(SEXP XsvtR, int nr, int j, double* col, int* idx) {
    SEXP Xsvt_SVT;
    SEXP svtLeaf;
    int  nnz = 0;

    Xsvt_SVT = GET_SLOT(XsvtR, SVT_SparseArray_svtSym);
    svtLeaf = VECTOR_ELT(Xsvt_SVT, j);

    if (svtLeaf != R_NilValue) {
        SEXP valsR = VECTOR_ELT(svtLeaf, 0);
        SEXP offsetsR = VECTOR_ELT(svtLeaf, 1);
        int  nvals = length(valsR);
        int  noffsets = length(offsetsR);

        if (nvals > 0) {
          Memcpy(col, REAL(valsR), (size_t) nvals);
          Memcpy(idx, INTEGER(offsetsR), (size_t) noffsets);
        } else { /* lacunar */
          for (int i=0; i < noffsets; i++)
            col[i] = 1;
        }
        nnz = noffsets;
    }

    return nnz;
}

typedef int (*FetchNzColFunDef)(SEXP, int, int, double*, int*);

FetchNzColFunDef
find_dim_and_fetchnzcolfun(SEXP XR, int** dim) {
  FetchNzColFunDef fetch_col;
  const char*      class = get_class_name(XR);

  if (!strcmp(class, "matrix")) {
    fetch_col = &fetchnzcol_dblmatrix;
    *dim = INTEGER(getAttrib(XR, R_DimSymbol));
  } else if (!strcmp(class, "dgCMatrix")) {
    fetch_col = &fetchnzcol_dgCMatrix;
    *dim = INTEGER(GET_SLOT(XR, Matrix_DimSym));
  } else if (!strcmp(class, "SVT_SparseMatrix")) {
    fetch_col = &fetchnzcol_dblSVT_SparseMatrix;
    *dim = INTEGER(GET_SLOT(XR, SVT_SparseArray_dimSym));
  } else
    error("input class %s cannot be handled yet.", class);

  return fetch_col;
}

/* assumes genesetsidxR is a list of ordered integer vectors, each containing
 * the indices of the genes in the corresponding gene set. */
SEXP
avg_score_genesets_R(SEXP XR, SEXP genesetsidxR,
                     SEXP anynaR, SEXP nauseR, SEXP minsizeR, SEXP verboseR) {
  int*     dimX;
  int      p, n;
  int      m = length(genesetsidxR);
  Rboolean anyna=asLogical(anynaR);
  int      nause=INTEGER(nauseR)[0]; /* everything=1, all.obs=2, na.rm=3 */
  int      minsize=INTEGER(minsizeR)[0];
  SEXP     esR;
  double*  es;
  int      wna=0;
  Rboolean abort=FALSE;
  Rboolean verbose=asLogical(verboseR);
  SEXP     pb=R_NilValue;
  int      nunprotect=0;
  double*  col;
  int*     off;
  FetchNzColFunDef fetch_col;

  fetch_col = find_dim_and_fetchnzcolfun(XR, &dimX);
  p = dimX[0]; /* number of rows/genes/features */
  n = dimX[1]; /* number of columns/samples/cells/spots */

  col = R_Calloc(p, double);
  off = R_Calloc(p, int);

  PROTECT(esR = allocMatrix(REALSXP, m, n)); nunprotect++;
  es = REAL(esR);

  if (verbose) {
    pb = PROTECT(cli_progress_bar(p, NULL)); nunprotect++;
    cli_progress_set_name(pb, "Calculating average scores");
  }

  for (int i=0; i < n; i++) {
    int nnz;

    if (verbose) { /* show progress */
      if (i % 100 == 0 && CLI_SHOULD_TICK)
        cli_progress_set(pb, i);
    }

    nnz = (*fetch_col)(XR, p, i, col, off);

    for (int j=0; j < m; j++) {
      SEXP     gsetidxR = VECTOR_ELT(genesetsidxR, j);
      int*     gsetidx;
      int      k = length(gsetidxR);
#ifdef LONG_VECTOR_SUPPORT
      R_xlen_t idx = (R_xlen_t) m * i + j;
#else
      int      idx = (size_t) m * i + j;
#endif
      double   sum = 0.0;
      int      idx_gs;
      int      idx_col;
      int      nnas = 0; /* number of NAs in the current gene set */

      gsetidx = INTEGER(gsetidxR); /* assume indices in gsetidx are ordered and 1-based */

      /* sum nonzero values in the column for the current gene set */
      if (nnz == p) { /* dense column */
        for (idx_gs=0; idx_gs < k; idx_gs++) {
          if (anyna && ISNA(col[gsetidx[idx_gs] - 1]) && nause == 2) { /* all.obs */
            abort=TRUE;
            break;
          } else if (anyna && ISNA(col[gsetidx[idx_gs] - 1]) && nause == 3) /* na.rm */
            nnas++;
          else
            sum = sum + col[gsetidx[idx_gs] - 1];
        }
      } else { /* sparse column */
        int prev_idx_gs = -1; /* previous gene set index, for checking order */

        idx_gs = 0;
        idx_col = 0;
        while (idx_gs < k && idx_col < nnz) {
          while (idx_col < nnz && gsetidx[idx_gs] - 1 > off[idx_col])
            idx_col++;
          if (idx_col < nnz && gsetidx[idx_gs] - 1 == off[idx_col]) {
            if (anyna && ISNA(col[idx_col]) && nause == 2) { /* all.obs */
              abort=TRUE;
              break;
            } else if (anyna && ISNA(col[idx_col]) && nause == 3) /* na.rm */
              nnas++;
            else
              sum = sum + col[idx_col];
          }

          if (prev_idx_gs >= 0 && gsetidx[idx_gs] <= prev_idx_gs)
            error("avg_score_genesets_R: gene set indices must be ordered and unique");
          prev_idx_gs = gsetidx[idx_gs];

          idx_gs++;
        }
      }
      es[idx] = NA_REAL;
      if (abort)
        break;
      else if (k - nnas >= minsize)
        es[idx] = sum / (double) (k - nnas);
      else
        wna = 1; /* warn about NAs in the output */
    }
  }

  R_Free(col);
  R_Free(off);

  if (anyna) {
    SEXP attr;

    if (nause == 2 && abort) {
      PROTECT(attr = allocVector(STRSXP, 1));
      SET_STRING_ELT(attr, 0, mkChar("abort"));
      Rf_setAttrib(esR, GSVA_attrNAsSym, attr);
      UNPROTECT(1); /* attr */
    } else if (nause == 3 && wna == 1) {
      PROTECT(attr = allocVector(STRSXP, 1));
      SET_STRING_ELT(attr, 0, mkChar("wna"));
      Rf_setAttrib(esR, GSVA_attrNAsSym, attr);
      UNPROTECT(1); /* attr */
    }
  }

  if (verbose)
    cli_progress_done(pb);

  UNPROTECT(nunprotect); /* esR pb */

  return(esR);
}

