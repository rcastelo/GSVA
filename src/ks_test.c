#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include <cli/progress.h>
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
gsva_sparse_rnd_walk(int* gsetidx, int k, int* decordstat, double* symrnkstat, int n,
                     double tau, int nzeros, double zerosymrnkstat,
                     double* walkstat, double* walkstatpos, double* walkstatneg) {
  int     nnz, gset_nnz, gset_nzeros;
  int*    gsetrnk;
  int*    gsetrnk_nnz;
  int*    gsetidx_nnz;
  double* stepcdfingeneset;
  int*    stepcdfoutgeneset;

  nnz = n - nzeros;

  gsetrnk = R_Calloc(k, int);
  gsetrnk_nnz = R_Calloc(k, int);
  gsetidx_nnz = R_Calloc(k, int);
  gset_nnz = gset_nzeros = 0;

  for (int i=0; i < k; i++) {
    gsetrnk[i] = decordstat[gsetidx[i]-1];
    if (gsetrnk[i] <= nnz) { /* assuming nonzeros expression values are positive */
      gsetrnk_nnz[gset_nnz] = gsetrnk[i];
      gsetidx_nnz[gset_nnz] = gsetidx[i];
      gset_nnz++;
    } else
      gset_nzeros++;
  }

  stepcdfingeneset = R_Calloc(n, double);  /* assuming zeros are set */
  stepcdfoutgeneset = R_Calloc(n, int);
  for (int i=0; i < n; i++)
    stepcdfoutgeneset[i] = 1;

  if (gset_nnz > 0) { /* if there are nonzero expression values in the gene set */
    for (int i=0; i < gset_nnz; i++) {
      /* convert 1-based gene indices to 0-based ! */
      if (tau == 1)
        stepcdfingeneset[gsetrnk_nnz[i]-1] = symrnkstat[gsetidx_nnz[i]-1];
      else
        stepcdfingeneset[gsetrnk_nnz[i]-1] = pow(symrnkstat[gsetidx_nnz[i]-1], tau);
      stepcdfoutgeneset[gsetrnk_nnz[i]-1] = 0;
    }
  }

  R_Free(gsetrnk);
  R_Free(gsetrnk_nnz);
  R_Free(gsetidx_nnz);

  if (gset_nzeros > 0) { /* if there are zero expression values in the gene set */
    int*   gsetrnk_zeros;
    double step = 0;
   
    if (gset_nzeros > 1)
      step = (double) (nzeros - 1) / (double) (gset_nzeros - 1);

    gsetrnk_zeros = R_Calloc(gset_nzeros, int);
    for (int i=0; i < gset_nzeros; i++)
      gsetrnk_zeros[i] = ((int) round(((double) i) * step + 1)) + nnz;

    for (int i=0; i < gset_nzeros; i++) {
      /* convert 1-based gene indices to 0-based ! */
      if (tau == 1)
        stepcdfingeneset[gsetrnk_zeros[i]-1] = zerosymrnkstat;
      else
        stepcdfingeneset[gsetrnk_zeros[i]-1] = pow(zerosymrnkstat, tau);
      stepcdfoutgeneset[gsetrnk_zeros[i]-1] = 0;
    }

    R_Free(gsetrnk_zeros);
  }

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

    for (int i=0; i < k_notna; i++) {
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
gsva_score_genesets_R(SEXP ranksR, SEXP genesetsidxR, SEXP intrnksR,
                      SEXP sparseR, SEXP maxdiffR, SEXP absrnkR, SEXP tauR,
                      SEXP anynaR, SEXP nauseR, SEXP minsizeR, SEXP verboseR) {
  int      p, n, nnas;
  int      m = length(genesetsidxR);
  Rboolean intrnks=asLogical(intrnksR);
  Rboolean sparse=asLogical(sparseR);
  Rboolean maxdiff=asLogical(maxdiffR);
  Rboolean absrnk=asLogical(absrnkR);
  double   tau=REAL(tauR)[0];
  Rboolean anyna=asLogical(anynaR);
  int      nause=INTEGER(nauseR)[0]; /* everything=1, all.obs=2, na.rm=3 */
  int      minsize=INTEGER(minsizeR)[0];
  SEXP     esR;
  double*  es;
  int      wna = 0;
  int      nzeros_col = 0;
  double   zerosymrnkstat = 0.0;
  int      invalidranks = 0;
  Rboolean abort = FALSE;
  Rboolean verbose = asLogical(verboseR);
  SEXP     pb = R_NilValue;
  int      nunprotect = 0;
  int*     decordstat_col;
  double*  symrnkstat_col;
  ranks_ctx_t* ctx;

  ctx = ranks_ctx_create(ranksR, intrnks, sparse);
  p = ctx->p; /* number of rows/genes */
  n = ctx->n; /* number of columns/samples */

  decordstat_col = (int*)R_alloc(p, sizeof(int));
  symrnkstat_col = (double*)R_alloc(p, sizeof(double));

  PROTECT(esR = allocMatrix(REALSXP, m, n)); nunprotect++;
  es = REAL(esR);

  if (verbose) {
    pb = PROTECT(cli_progress_bar(n, NULL)); nunprotect++; /* one tick per column */
    cli_progress_set_name(pb, "Calculating GSVA scores");
  }

  for (int i=0; i < n; i++) {
    if (verbose) { /* show progress */
      if (i % 100 == 0 && CLI_SHOULD_TICK)
        cli_progress_set(pb, i);
    }

    if (anyna)
      invalidranks = ranks2stats_nas(ctx, i, decordstat_col, symrnkstat_col,
                                     &nzeros_col, &zerosymrnkstat, &nnas);
    else
      invalidranks = ranks2stats(ctx, i, decordstat_col, symrnkstat_col,
                                 &nzeros_col, &zerosymrnkstat);

    /* ranks out of range would lead to reading and writing out of the
     * bounds of the arrays used in the random walk */
    if (invalidranks) {
      if (verbose)
        cli_progress_done(pb);
      UNPROTECT(nunprotect); /* esR pb */
      if (invalidranks == 2)
        error("GSVA ranks have missing values, but the GSVA parameters "
              "indicate that the input data has none, e.g., because "
              "'checkNA=\"no\"' was used in 'gsvaParam()'.");
      else
        error("GSVA ranks are out of range for the %d rows of the input "
              "data. This may happen when rows are removed after "
              "calculating the ranks with 'gsvaColRanks()'; remove them "
              "before calling 'gsvaRowNorm()' instead.", p);
    }

    for (int j=0; j < m; j++) {
      SEXP     gsetidxR = VECTOR_ELT(genesetsidxR, j);
      int*     gsetidx;
      int      k = length(gsetidxR);
#ifdef LONG_VECTOR_SUPPORT
      R_xlen_t idx = (R_xlen_t) m * i + j;
#else
      int      idx = (size_t) m * i + j;
#endif
      double  walkstatpos, walkstatneg;

      walkstatpos = walkstatneg = NA_REAL;
      gsetidx = INTEGER(gsetidxR);
      if (anyna) {
        gsva_rnd_walk_nas(gsetidx, k, decordstat_col, symrnkstat_col, p, tau,
                          nause, minsize, NULL, &walkstatpos, &walkstatneg, &wna);
      } else {
        if (sparse && nzeros_col > 0)
          gsva_sparse_rnd_walk(gsetidx, k, decordstat_col, symrnkstat_col, p,
                               tau, nzeros_col, zerosymrnkstat, NULL,
                               &walkstatpos, &walkstatneg);
        else
          gsva_rnd_walk(gsetidx, k, decordstat_col, symrnkstat_col, p, tau,
                        NULL, &walkstatpos, &walkstatneg);
      }

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
    if (abort)
      break;
  }

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

/* figures out whether 'obj' is a matrix, dgCMatrix or SVT_SparseMatrix
 * and returns the class name */
const char*
get_class_name(SEXP obj) {
    if (isMatrix(obj))
        return "matrix";

    SEXP classR = getAttrib(obj, R_ClassSymbol);
    return CHAR(STRING_ELT(classR, 0));
}
