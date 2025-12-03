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
            SVT_SparseArray_svtSym;

/* to add attributes to R objects from C code */
static SEXP
installAttrib(SEXP, SEXP, SEXP);

static SEXP
installAttrib(SEXP vec, SEXP name, SEXP val)
{
  SEXP s, t;

  if (TYPEOF(vec) == CHARSXP)
    error("cannot set attribute on a CHARSXP");
  PROTECT(vec);
  PROTECT(name);
  PROTECT(val);
  for (s = ATTRIB(vec); s != R_NilValue; s = CDR(s)) {
    if (TAG(s) == name) {
      SETCAR(s, val);
      UNPROTECT(3);
      return val;
    }
  }
  s = Rf_allocList(1);
  SETCAR(s, val);
  SET_TAG(s, name);
  if (ATTRIB(vec) == R_NilValue)
    SET_ATTRIB(vec, s);
  else {
    t = nthcdr(ATTRIB(vec), length(ATTRIB(vec)) - 1);
    SETCDR(t, s);
  }
  UNPROTECT(3);
  return val;
}

void
gsva_rnd_walk(int* gsetidx, int k, int* decordstat, double* symrnkstat, int n,
              double tau, double* walkstat, double* walkstatpos,
              double* walkstatneg) {
  int*    gsetrnk;
  double* stepcdfingeneset;
  int*    stepcdfoutgeneset;

  gsetrnk = R_Calloc(k, int);
  for (int i=0; i < k; i++)
    gsetrnk[i] = decordstat[gsetidx[i]-1];

  stepcdfingeneset = R_Calloc(n, double);  /* assuming zeroes are set */
  stepcdfoutgeneset = R_Calloc(n, int);
  for (int i=0; i < n; i++)
    stepcdfoutgeneset[i] = 1;

  for (int i=0; i < k; i++) {
    /* convert 1-based gene indices to 0-based ! */
    if (tau == 1)
      stepcdfingeneset[gsetrnk[i]-1] = symrnkstat[gsetidx[i]-1];
    else
      stepcdfingeneset[gsetrnk[i]-1] = pow(symrnkstat[gsetidx[i]-1], tau);
    stepcdfoutgeneset[gsetrnk[i]-1] = 0;
  }
  R_Free(gsetrnk);

  for (int i=1; i < n; i++) {
    stepcdfingeneset[i] = stepcdfingeneset[i-1] + stepcdfingeneset[i];
    stepcdfoutgeneset[i] = stepcdfoutgeneset[i-1] + stepcdfoutgeneset[i];
  }

  *walkstatpos = *walkstatneg = NA_REAL;
  if (stepcdfingeneset[n-1] > 0 && stepcdfoutgeneset[n-1] > 0) {
    *walkstatpos = *walkstatneg = 0;
    for (int i=0; i < n; i++) {
      double wlkstat = 0;

      if (walkstat != NULL)
        wlkstat = walkstat[i] = ((double) stepcdfingeneset[i]) /
                                ((double) stepcdfingeneset[n-1]) -
                                ((double) stepcdfoutgeneset[i]) /
                                ((double) stepcdfoutgeneset[n-1]);
      else {
        wlkstat = ((double) stepcdfingeneset[i]) / ((double) stepcdfingeneset[n-1]) -
                  ((double) stepcdfoutgeneset[i]) / ((double) stepcdfoutgeneset[n-1]);
      }

      if (wlkstat > *walkstatpos)
        *walkstatpos = wlkstat;
      if (wlkstat < *walkstatneg)
        *walkstatneg = wlkstat;
    }
  }

  R_Free(stepcdfoutgeneset);
  R_Free(stepcdfingeneset);
}

void
gsva_rnd_walk_nas(int* gsetidx, int k, int* decordstat, double* symrnkstat, int n,
                  double tau, int na_use, int minsize, double* walkstat,
                  double* walkstatpos, double* walkstatneg, int* wna) {
  int*    gsetidx_wonas;
  int*    gsetrnk;
  double* stepcdfingeneset;
  int*    stepcdfoutgeneset;
  int     k_notna = 0;

  gsetidx_wonas = R_Calloc(k, int);
  gsetrnk = R_Calloc(k, int);

  for (int i=0; i < k; i++) {
    if (decordstat[gsetidx[i]-1] != NA_INTEGER) { /* na.rm skips NAs */
      gsetidx_wonas[k_notna] = gsetidx[i];
      gsetrnk[k_notna] = decordstat[gsetidx[i]-1];
      k_notna++;
    } else {
      if (na_use < 3) /* everything or all.obs */
        return;
    }
  }

  *walkstatpos = *walkstatneg = NA_REAL;
  if (k_notna >= minsize) { /* na.rm */
    k = k_notna;

    stepcdfingeneset = R_Calloc(n, double);  /* assuming zeroes are set */
    stepcdfoutgeneset = R_Calloc(n, int);
    for (int i=0; i < n; i++)
      stepcdfoutgeneset[i] = 1;

    for (int i=0; i < k; i++) {
      /* convert 1-based gene indices to 0-based ! */
      if (tau == 1)
        stepcdfingeneset[gsetrnk[i]-1] = symrnkstat[gsetidx_wonas[i]-1];
      else
        stepcdfingeneset[gsetrnk[i]-1] = pow(symrnkstat[gsetidx_wonas[i]-1], tau);
      stepcdfoutgeneset[gsetrnk[i]-1] = 0;
    }

    for (int i=1; i < n; i++) {
      stepcdfingeneset[i] = stepcdfingeneset[i-1] + stepcdfingeneset[i];
      stepcdfoutgeneset[i] = stepcdfoutgeneset[i-1] + stepcdfoutgeneset[i];
    }

    if (stepcdfingeneset[n-1] > 0 && stepcdfoutgeneset[n-1] > 0) {
      *walkstatpos = *walkstatneg = 0;
      for (int i=0; i < n; i++) {
        double wlkstat = 0;

        if (walkstat != NULL)
          wlkstat = walkstat[i] = ((double) stepcdfingeneset[i]) /
                                  ((double) stepcdfingeneset[n-1]) -
                                  ((double) stepcdfoutgeneset[i]) /
                                  ((double) stepcdfoutgeneset[n-1]);
        else {
          wlkstat = ((double) stepcdfingeneset[i]) /
                    ((double) stepcdfingeneset[n-1]) -
                    ((double) stepcdfoutgeneset[i]) /
                    ((double) stepcdfoutgeneset[n-1]);
        }

        if (wlkstat > *walkstatpos)
          *walkstatpos = wlkstat;
        if (wlkstat < *walkstatneg)
          *walkstatneg = wlkstat;
      }
    }

    R_Free(stepcdfoutgeneset);
    R_Free(stepcdfingeneset);

  } else
    *wna = 1;

  R_Free(gsetrnk);
  R_Free(gsetidx_wonas);
}

SEXP
old_gsva_score_genesets_R(SEXP genesetsidxR, SEXP decordstatR, SEXP symrnkstatR,
                          SEXP maxdiffR, SEXP absrnkR, SEXP tauR, SEXP anynaR,
                          SEXP nauseR, SEXP minsizeR, SEXP verboseR) {
  SEXP     dimInput;
  int      m = length(genesetsidxR);
  int      p, n;
  Rboolean maxdiff=asLogical(maxdiffR);
  Rboolean absrnk=asLogical(absrnkR);
  double   tau=REAL(tauR)[0];
  Rboolean anyna=asLogical(anynaR);
  int      nause=INTEGER(nauseR)[0]; /* everything=1, all.obs=2, na.rm=3 */
  int      minsize=INTEGER(minsizeR)[0];
  int*     decordstat;
  double*  symrnkstat;
  SEXP     esR;
  double*  es;
  int      wna=0;
  Rboolean abort=FALSE;
  Rboolean verbose=asLogical(verboseR);
  SEXP     pb=R_NilValue;
  int      nunprotect=0;

  dimInput = getAttrib(decordstatR, R_DimSymbol);
  p = INTEGER(dimInput)[0]; /* number of genes/features */
  n = INTEGER(dimInput)[1]; /* number of samples/cells */

  decordstat = INTEGER(decordstatR);
  symrnkstat = REAL(symrnkstatR);

  PROTECT(esR = allocMatrix(REALSXP, m, n)); nunprotect++;
  es = REAL(esR);

  if (verbose) {
    pb = PROTECT(cli_progress_bar(p, NULL)); nunprotect++;
    cli_progress_set_name(pb, "Calculating GSVA scores");
  }

  for (int i=0; i < n; i++) {
    int*     decordstat_col = decordstat + i * p;
    double*  symrnkstat_col = symrnkstat + i * p;

    if (verbose) { /* show progress */
      if (i % 100 == 0 && CLI_SHOULD_TICK)
        cli_progress_set(pb, i);
    }

    for (int j=0; j < m; j++) {
      SEXP    gsetidxR=VECTOR_ELT(genesetsidxR, j);
      int*    gsetidx;
      int     k = length(gsetidxR);
#ifdef LONG_VECTOR_SUPPORT
      R_xlen_t idx = m * i + j;
#else
      int idx = m * i + j;
#endif
      double  walkstatpos, walkstatneg;

      walkstatpos = walkstatneg = NA_REAL;
      gsetidx = INTEGER(gsetidxR);
      if (anyna)
        gsva_rnd_walk_nas(gsetidx, k, decordstat_col, symrnkstat_col, p, tau,
                          nause, minsize, NULL, &walkstatpos, &walkstatneg,
                          &wna);
      else
        gsva_rnd_walk(gsetidx, k, decordstat_col, symrnkstat_col, p, tau,
                      NULL, &walkstatpos, &walkstatneg);

      es[idx] = NA_REAL;
      if (!anyna || (!ISNA(walkstatpos) && !ISNA(walkstatneg))) {
	      if (maxdiff) {
		      es[idx] = walkstatpos + walkstatneg;
          if (absrnk)
            es[idx] = walkstatpos - walkstatneg;
	      } else {
		        es[idx] = (walkstatpos > fabs(walkstatneg)) ? walkstatpos : walkstatneg;
	      }
      } else {
        if (anyna && (ISNA(walkstatpos) || ISNA(walkstatneg)) && nause == 2) { /* all.obs */
          abort=TRUE;
          break;
        }
      }
    }
  }

  if (anyna) {
    SEXP class;

    if (nause == 2 && abort) {
      PROTECT(class = allocVector(STRSXP, 1));
      SET_STRING_ELT(class, 0, mkChar("abort"));
      installAttrib(esR, R_ClassSymbol, class);
      UNPROTECT(1); /* class */
    } else if (nause == 3 && wna == 1) {
      PROTECT(class = allocVector(STRSXP, 1));
      SET_STRING_ELT(class, 0, mkChar("wna"));
      installAttrib(esR, R_ClassSymbol, class);
      UNPROTECT(1); /* class */
    }
  }

  if (verbose)
    cli_progress_done(pb);

  UNPROTECT(nunprotect); /* esR pb */

  return(esR);
}

/* fetch column from a dense matrix XR
 * nr - number of rows
 * j - 0-based column to fetch
 * col - array where to store the column, assuming is initialized to zeroes
 * returned value - number of rows
 */
int
fetch_col_matrix(SEXP XR, int nr, int j, int* col) {
  int* X=INTEGER(XR);

  Memcpy(col, X+nr*j, (size_t) nr);

  return nr;
}

/* fetch column from a sparse dgCMatrix XR
 * nr - number of rows
 * j - 0-based column to fetch
 * col - array where to store the column, assuming is initialized to zeroes
 * returned value - number of nonzero values
 */
int
fetch_col_dgCMatrix(SEXP XCspR, int nr, int j, int* col) {
  int*    XCsp_i;
  int*    XCsp_p;
  double* XCsp_x;

  XCsp_i = INTEGER(GET_SLOT(XCspR, Matrix_iSym));
  XCsp_p = INTEGER(GET_SLOT(XCspR, Matrix_pSym));
  XCsp_x = REAL(GET_SLOT(XCspR, Matrix_xSym));

  /* put the sparse column into a dense vector */
  for (int i=XCsp_p[j]; i < XCsp_p[j+1]; i++)
    col[XCsp_i[i]] = (int) XCsp_x[i];

  return XCsp_p[j+1]-XCsp_p[j];
}

/* fetch column from a sparse SVT_SparseMatrix XR
 * nr - number of rows
 * j - 0-based column to fetch
 * col - array where to store the column, assuming is initialized to zeroes
 * returned value - number of nonzero values
 */
int
fetch_col_SVT_SparseMatrix(SEXP XsvtR, int nr, int j, int* col) {
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

typedef int (*FetchColFunDef)(SEXP, int, int, int*);

FetchColFunDef
find_dim_and_fetchcolfun(SEXP XR, int** dim) {
  FetchColFunDef fetch_col;
  SEXP           classR = eval(lang2(install("class"), XR), R_GlobalEnv);
  const char*    class = CHAR(STRING_ELT(classR, 0));

  if (!strcmp(class, "matrix")) {
    fetch_col = &fetch_col_matrix;
    *dim = INTEGER(getAttrib(XR, R_DimSymbol));
  } else if (!strcmp(class, "dgCMatrix")) {
    fetch_col = &fetch_col_dgCMatrix;
    *dim = INTEGER(GET_SLOT(XR, Matrix_DimSym));
  } else if (!strcmp(class, "SVT_SparseMatrix")) {
    fetch_col = &fetch_col_SVT_SparseMatrix;
    *dim = INTEGER(GET_SLOT(XR, SVT_SparseArray_dimSym));
  } else
    error("input class %s cannot be handled yet\n", class);

  return fetch_col;
}

void
ranks2stats(SEXP ranksR, int p, int n, int j, Rboolean sparse,
            FetchColFunDef fetch_col,
            int* decordstat_col, double* symrnkstat_col);

void
ranks2stats_nas(SEXP ranksR, int p, int n, int j, Rboolean sparse,
                FetchColFunDef fetch_col,
                int* decordstat_col, double* symrnkstat_col);

SEXP
gsva_score_genesets_R(SEXP ranksR, SEXP genesetsidxR, SEXP sparseR,
                      SEXP maxdiffR, SEXP absrnkR, SEXP tauR, SEXP anynaR,
                      SEXP nauseR, SEXP minsizeR, SEXP verboseR) {
  int*     dimranks;
  int      p, n;
  int      m = length(genesetsidxR);
  Rboolean sparse=asLogical(sparseR);
  Rboolean maxdiff=asLogical(maxdiffR);
  Rboolean absrnk=asLogical(absrnkR);
  double   tau=REAL(tauR)[0];
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
  int*     decordstat_col;
  double*  symrnkstat_col;
  FetchColFunDef fetch_col;

  fetch_col = find_dim_and_fetchcolfun(ranksR, &dimranks);
  p = dimranks[0]; /* number of rows/genes/features */
  n = dimranks[1]; /* number of columns/samples/cells/spots */

  decordstat_col = R_Calloc(p, int);
  symrnkstat_col = R_Calloc(p, double);

  PROTECT(esR = allocMatrix(REALSXP, m, n)); nunprotect++;
  es = REAL(esR);

  if (verbose) {
    pb = PROTECT(cli_progress_bar(p, NULL)); nunprotect++;
    cli_progress_set_name(pb, "Calculating GSVA scores");
  }

  for (int i=0; i < n; i++) {
    if (verbose) { /* show progress */
      if (i % 100 == 0 && CLI_SHOULD_TICK)
        cli_progress_set(pb, i);
    }

    if (anyna)
      ranks2stats_nas(ranksR, p, n, i, sparse, fetch_col, decordstat_col, symrnkstat_col);
    else
      ranks2stats(ranksR, p, n, i, sparse, fetch_col, decordstat_col, symrnkstat_col);

    for (int j=0; j < m; j++) {
      SEXP    gsetidxR=VECTOR_ELT(genesetsidxR, j);
      int*    gsetidx;
      int     k = length(gsetidxR);
#ifdef LONG_VECTOR_SUPPORT
      R_xlen_t idx = m * i + j;
#else
      int idx = m * i + j;
#endif
      double  walkstatpos, walkstatneg;

      walkstatpos = walkstatneg = NA_REAL;
      gsetidx = INTEGER(gsetidxR);
      if (anyna)
        gsva_rnd_walk_nas(gsetidx, k, decordstat_col, symrnkstat_col, p, tau,
                          nause, minsize, NULL, &walkstatpos, &walkstatneg,
                          &wna);
      else
        gsva_rnd_walk(gsetidx, k, decordstat_col, symrnkstat_col, p, tau,
                      NULL, &walkstatpos, &walkstatneg);

      es[idx] = NA_REAL;
      if (!anyna || (!ISNA(walkstatpos) && !ISNA(walkstatneg))) {
	      if (maxdiff) {
		      es[idx] = walkstatpos + walkstatneg;
          if (absrnk)
            es[idx] = walkstatpos - walkstatneg;
	      } else {
		        es[idx] = (walkstatpos > fabs(walkstatneg)) ? walkstatpos : walkstatneg;
	      }
      } else {
        if (anyna && (ISNA(walkstatpos) || ISNA(walkstatneg)) && nause == 2) { /* all.obs */
          abort=TRUE;
          break;
        }
      }
    }
  }

  R_Free(decordstat_col);
  R_Free(symrnkstat_col);

  if (anyna) {
    SEXP class;

    if (nause == 2 && abort) {
      PROTECT(class = allocVector(STRSXP, 1));
      SET_STRING_ELT(class, 0, mkChar("abort"));
      installAttrib(esR, R_ClassSymbol, class);
      UNPROTECT(1); /* class */
    } else if (nause == 3 && wna == 1) {
      PROTECT(class = allocVector(STRSXP, 1));
      SET_STRING_ELT(class, 0, mkChar("wna"));
      installAttrib(esR, R_ClassSymbol, class);
      UNPROTECT(1); /* class */
    }
  }

  if (verbose)
    cli_progress_done(pb);

  UNPROTECT(nunprotect); /* esR pb */

  return(esR);
}

/* j is a 0-based column index on ranksR */
void
ranks2stats(SEXP ranksR, int p, int n, int j, Rboolean sparse,
            FetchColFunDef fetch_col,
            int* decordstat_col, double* symrnkstat_col) {
  int* r = R_Calloc(p, int);       /* assume 0s are set */
  int* r_dense = R_Calloc(p, int); /* assume 0s are set */
  int  nnz, nzs;

  nnz = (*fetch_col)(ranksR, p, j, r);
  nzs = p - nnz;

  if (nzs > 0) { /* if ranks have zeroes, then input is a sparse matrix */
    int  k = 1;

    for (int i=0; i < p; i++) {
      if (r[i] == 0)       /* sparse ranks into dense ranks */
        r_dense[i] = k++;
      else
        r_dense[i] = r[i] + nzs;
    }
  } else         /* input is a dense matrix */
    for (int i=0; i < p; i++)
      r_dense[i] = r[i];

  /* dense ranks into decreasing order statistics */
  for (int i=0; i < p; i++)
    decordstat_col[i] = p - r_dense[i] + 1;

  if (nzs > 0 && sparse) {
    for (int i=0; i < p; i++) {
      double nnz1div2 = ((double) (nnz+1)) / 2.0;  /* nnz is in fact max(r)   */
      if (r[i] == 0)                               /* in sparse regime zeroes */
        symrnkstat_col[i] = fabs(nnz1div2 - 1.0);  /* get same sym rank stat  */
      else                                         /* nonzero ranks shift one */
        symrnkstat_col[i] = fabs(nnz1div2 - (double) (r[i] + 1));
    }
  } else {
    for (int i=0; i < p; i++)
      symrnkstat_col[i] = fabs(((double) p) / 2.0 - ((double) r_dense[i]));
  }

  R_Free(r_dense);
  R_Free(r);
}

/* j is a 0-based column index on ranksR
 * regular reminder that ISNA() only applies to numeric values of
 * type double and missingess of integer values should be tested
 * by equality to NA_INTEGER, see
 * https://cran.r-project.org/doc/manuals/r-devel/R-exts.html#Missing-and-special-values-1
 */
void
ranks2stats_nas(SEXP ranksR, int p, int n, int j, Rboolean sparse,
                FetchColFunDef fetch_col,
                int* decordstat_col, double* symrnkstat_col) {
  int* r = R_Calloc(p, int);       /* assume 0s are set */
  int* r_dense = R_Calloc(p, int); /* assume 0s are set */
  int  nnz, nzs;
  int  nnas;

  nnz = (*fetch_col)(ranksR, p, j, r);
  nnas = 0;
  for (int i=0; i < p; i++)
    if (r[i] == NA_INTEGER)
      nnas++;

  nzs = p - nnz;

  if (nzs > 0) { /* if ranks have zeroes, then input is a sparse matrix */
    int  k = 1;

    for (int i=0; i < p; i++) {
      if (r[i] != NA_INTEGER) {
        if (r[i] == 0)       /* sparse ranks into dense ranks */
          r_dense[i] = k++;
        else
          r_dense[i] = r[i] + nzs;
      } else
        r_dense[i] = NA_INTEGER;
    }
  } else         /* input is a dense matrix */
    for (int i=0; i < p; i++)
      r_dense[i] = r[i];

  /* dense ranks into decreasing order statistics */
  for (int i=0; i < p; i++)
    decordstat_col[i] = r_dense[i] == NA_INTEGER ? NA_INTEGER : p - nnas - r_dense[i] + 1;

  if (nzs > 0 && sparse) {
    for (int i=0; i < p; i++) {
      double nnz1div2 = ((double) (nnz-nnas+1)) / 2.0; /* nnz is in fact max(r) */
      if (r[i] != NA_INTEGER) {
        if (r[i] == 0)                               /* in sparse regime zeroes */
          symrnkstat_col[i] = fabs(nnz1div2 - 1.0);  /* get same sym rank stat  */
        else                                         /* nonzero ranks shift one */
          symrnkstat_col[i] = fabs(nnz1div2 - (double) (r[i] + 1));
      } else
          symrnkstat_col[i] = NA_REAL;
    }
  } else {
    for (int i=0; i < p; i++)
      if (r_dense[i] != NA_INTEGER)
        symrnkstat_col[i] = fabs(((double) (p - nnas)) / 2.0 - ((double) r_dense[i]));
      else
        symrnkstat_col[i] = NA_REAL;
  }

  R_Free(r_dense);
  R_Free(r);
}

/* R wrapper for testing purposes */
SEXP
ranks2stats_R(SEXP ranksR, SEXP jR, SEXP sparseR, SEXP anynaR) {
  int            j = INTEGER(jR)[0];
  Rboolean       sparse = asLogical(sparseR);
  Rboolean       anyna = asLogical(anynaR);
  int            p, n;
  int*           dimranks;
  FetchColFunDef fetch_col;
  SEXP           dosR, srsR, ansR;
  int*           dos;
  double*        srs;

  fetch_col = find_dim_and_fetchcolfun(ranksR, &dimranks);
  p = dimranks[0]; /* number of rows/genes/features */
  n = dimranks[1]; /* number of columns/samples/cells/spots */

  PROTECT(dosR = allocVector(INTSXP, p));
  PROTECT(srsR = allocVector(REALSXP, p));

  dos = INTEGER(dosR);
  srs = REAL(srsR);

  if (anyna)
    ranks2stats_nas(ranksR, p, n, j-1, sparse, fetch_col, dos, srs);
  else
    ranks2stats(ranksR, p, n, j-1, sparse, fetch_col, dos, srs);

  PROTECT(ansR = allocVector(VECSXP, 2));

  SET_VECTOR_ELT(ansR, 0, dosR);
  SET_VECTOR_ELT(ansR, 1, srsR);

  UNPROTECT(3); /* dosR srsR ansR */

  return(ansR);
}
