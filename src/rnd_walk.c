#include <stdlib.h>
#include <math.h>
#include <R_ext/Utils.h>
#include "rnd_walk.h"

/* memory is allocated with R_alloc(), which is released by R when the .Call()
 * returns, also on error()
 */
rndwalk_ctx_t*
rndwalk_ctx_create(SEXP genesetsidxR) {
  int            m = length(genesetsidxR);
  int            maxk = 0;
  rndwalk_ctx_t* ctx;

  for (int j=0; j < m; j++) {
    int k = length(VECTOR_ELT(genesetsidxR, j));
    if (k > maxk)
      maxk = k;
  }

  ctx = (rndwalk_ctx_t*) R_alloc(1, sizeof(rndwalk_ctx_t));
  ctx->maxk = maxk;
  ctx->ranks = (int*) R_alloc(maxk, sizeof(int));
  ctx->idx = (int*) R_alloc(maxk, sizeof(int));
  ctx->weights = (double*) R_alloc(maxk, sizeof(double));

  return ctx;
}

/* maximum and minimum deviations of a random walk of length 'p' that goes up
 * at the 'k' positions in 'ranks', with the 'weights' that sum up to
 * 'total_pos_sum', and goes down by 1/(p-k) at every other position. 'ranks'
 * should hold 0-based, distinct, positions sorted in increasing order, and
 * 'idx' should give, for the i-th of them, the index in 'weights' of its
 * weight.
 *
 * the deviations are computed without walking through the 'p' positions: the
 * walk goes up only at the positions in 'ranks' and goes down at every other
 * position, so its maximum is attained right after a step up and its minimum
 * right before one. at the i-th position in 'ranks' (0-based) there are
 * (ranks[i] - i) steps down before it, which gives the value of the walk
 * right before (current_valley) and right after (current_peak) its step up.
 * walking through the 'k' positions instead of the 'p' ones reduces the
 * computational cost from O(n) to O(k). the caller should check beforehand
 * that the walk is defined, i.e., 0 < k < p and total_pos_sum > 0.
 */
static void
rnd_walk_extremes(const int* ranks, const int* idx, const double* weights,
                  int k, int p, double total_pos_sum, double* walkstatpos,
                  double* walkstatneg) {
  double running_pos_sum = 0.0;
  double max_peak = 0.0;
  double min_valley = 0.0;
  double neg_step = 1.0 / ((double) (p - k));

  for (int i=0; i < k; i++) {
    double neg_penalty = ((double) (ranks[i] - i)) * neg_step;
    double current_valley, current_peak;

    current_valley = running_pos_sum / total_pos_sum - neg_penalty;
    if (current_valley < min_valley)
      min_valley = current_valley;

    running_pos_sum += weights[idx[i]];

    current_peak = running_pos_sum / total_pos_sum - neg_penalty;
    if (current_peak > max_peak)
      max_peak = current_peak;
  }

  *walkstatpos = max_peak;
  *walkstatneg = min_valley;
}

/* the random walk is undefined, and its maximum deviations are set to NA,
 * when its in-gene-set or its out-gene-set steps sum up to zero. the latter
 * happens when the gene set contains all (nonmissing) genes and the former
 * when all genes in the gene set have a zero symmetric rank statistic, e.g.,
 * a single gene at the middle rank p/2 of a dense column with p even, or only
 * zero-expressed genes in a sparse column with a single nonzero value.
 *
 * the positions of the genes of the gene set are sorted with R_qsort_int_I(),
 * which also permutes the indices of the genes to read their weights in that
 * order, and then rnd_walk_extremes() computes the maximum and minimum
 * deviations of the walk. this reduces the original O(p) computational cost to
 * O(k log k) time. note that gene indices in 'gsetidx' are assumed to be
 * unique, which is guaranteed by the R code that maps gene sets to features.
 *
 * returned value: 1 if the random walk is undefined, 0 otherwise
 */
int
gsva_rnd_walk(rndwalk_ctx_t* ctx, int* gsetidx, int k, int* decordstat,
              double* symrnkstat, int p, double tau, double* walkstatpos,
              double* walkstatneg) {
  int*    ranks = ctx->ranks;
  int*    idx = ctx->idx;
  double* weights = ctx->weights;
  double  total_pos_sum = 0.0;

  *walkstatpos = *walkstatneg = NA_REAL;
  if (k == 0 || k == p) /* no in-gene-set or no out-gene-set steps */
    return 1;

  for (int i=0; i < k; i++) {
    /* convert 1-based gene indices and positions to 0-based ! */
    int gene_idx = gsetidx[i] - 1;

    ranks[i] = decordstat[gene_idx] - 1;
    idx[i] = i; /* reset on every call: R_qsort_int_I() permutes 'idx' in
                 * place, and it must start as 0..k-1 to pair each sorted
                 * position in 'ranks' with its weight in 'weights' */
    weights[i] = tau == 1 ? symrnkstat[gene_idx] :
                           pow(symrnkstat[gene_idx], tau);
    total_pos_sum += weights[i];
  }

  if (total_pos_sum == 0.0) /* in-gene-set steps sum up to zero */
    return 1;

  /* R_qsort_int_I() takes 1-based inclusive bounds, so 1, k sorts
   * ranks[0..k-1] */
  R_qsort_int_I(ranks, idx, 1, k);

  rnd_walk_extremes(ranks, idx, weights, k, p, total_pos_sum, walkstatpos,
                    walkstatneg);

  return 0;
}

/* returned value: 1 if the random walk is undefined, 0 otherwise
 * see gsva_rnd_walk()
 */
int
gsva_sparse_rnd_walk(int* gsetidx, int k, int* decordstat, double* symrnkstat,
                     int n, double tau, int nzeros, double zerosymrnkstat,
                     double* walkstatpos, double* walkstatneg) {
  int     nnz, gset_nnz, gset_nzeros;
  int*    gsetrnk;
  int*    gsetrnk_nnz;
  int*    gsetidx_nnz;
  double* stepcdfingeneset;
  int*    stepcdfoutgeneset;
  int     undefined;

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
  undefined = !(stepcdfingeneset[n-1] > 0 && stepcdfoutgeneset[n-1] > 0);
  if (!undefined) {
    *walkstatpos = *walkstatneg = 0;
    for (int i=0; i < n; i++) {
      double wlkstat = ((double) stepcdfingeneset[i]) /
                       ((double) stepcdfingeneset[n-1]) -
                       ((double) stepcdfoutgeneset[i]) /
                       ((double) stepcdfoutgeneset[n-1]);

      if (wlkstat > *walkstatpos)
        *walkstatpos = wlkstat;
      if (wlkstat < *walkstatneg)
        *walkstatneg = wlkstat;
    }
  }

  R_Free(stepcdfoutgeneset);
  R_Free(stepcdfingeneset);

  return undefined;
}

/* 'n' is the number of nonmissing genes in the column, i.e., the length of
 * the random walk, since nonmissing genes have decreasing order statistics
 * between 1 and 'n', and missing genes are left out of the walk.
 * returned value: 1 if the random walk is undefined, 0 otherwise, see
 * gsva_rnd_walk(); a gene set with missing values under na_use < 3 or with
 * less than 'minsize' nonmissing genes does not lead to an undefined walk
 */
int
gsva_rnd_walk_nas(int* gsetidx, int k, int* decordstat, double* symrnkstat,
                  int n, double tau, int na_use, int minsize,
                  double* walkstatpos, double* walkstatneg, int* wna) {
  int*    gsetidx_wonas;
  int*    gsetrnk;
  double* stepcdfingeneset;
  int*    stepcdfoutgeneset;
  int     k_notna = 0;
  int     undefined = 0;

  gsetidx_wonas = R_Calloc(k, int);
  gsetrnk = R_Calloc(k, int);

  for (int i=0; i < k; i++) {
    if (decordstat[gsetidx[i]-1] != NA_INTEGER) { /* na.rm skips NAs */
      gsetidx_wonas[k_notna] = gsetidx[i];
      gsetrnk[k_notna] = decordstat[gsetidx[i]-1];
      k_notna++;
    } else {
      if (na_use < 3) { /* everything or all.obs */
        R_Free(gsetrnk);
        R_Free(gsetidx_wonas);
        return 0;
      }
    }
  }

  *walkstatpos = *walkstatneg = NA_REAL;
  if (k_notna >= minsize) { /* na.rm */
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

    undefined = !(stepcdfingeneset[n-1] > 0 && stepcdfoutgeneset[n-1] > 0);
    if (!undefined) {
      *walkstatpos = *walkstatneg = 0;
      for (int i=0; i < n; i++) {
        double wlkstat = ((double) stepcdfingeneset[i]) /
                         ((double) stepcdfingeneset[n-1]) -
                         ((double) stepcdfoutgeneset[i]) /
                         ((double) stepcdfoutgeneset[n-1]);

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

  return undefined;
}
