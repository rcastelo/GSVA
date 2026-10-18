#include <stdlib.h>
#include <math.h>
#include <R.h>
#include "rnd_walk.h"

static int compare_gene_steps(const void *a, const void *b) {
    gene_step_t *stepA = (gene_step_t *)a;
    gene_step_t *stepB = (gene_step_t *)b;
    return (stepA->rank - stepB->rank);
}

rndwalk_ctx_t* rndwalk_ctx_create(SEXP genesetsidxR) {
    int            m = length(genesetsidxR);
    int            maxk = 0;
    rndwalk_ctx_t* ctx;

    for (int j = 0; j < m; j++) {
        R_xlen_t len = length(VECTOR_ELT(genesetsidxR, j));
        if (len > maxk) maxk = (int)len;
    }

    ctx = (rndwalk_ctx_t*)R_alloc(1, sizeof(rndwalk_ctx_t));
    ctx->maxk  = maxk;
    ctx->steps = (gene_step_t*)R_alloc(maxk, sizeof(gene_step_t));
    return ctx;
}

void gsva_rnd_walk(rndwalk_ctx_t* ctx, int* gsetidx, int k, int* decordstat,
                   double* symrnkstat, int n, double tau, double* walkstatpos,
                   double* walkstatneg) {
    if (k == 0 || k == n) {
        *walkstatpos = 0.0;
        *walkstatneg = 0.0;
        return;
    }

    double total_pos_sum = 0.0;

    for (int i = 0; i < k; i++) {
        int gene_idx = gsetidx[i];
        ctx->steps[i].rank = decordstat[gene_idx - 1] - 1;
        ctx->steps[i].score = pow(fabs(symrnkstat[gene_idx - 1]), tau);
        total_pos_sum += ctx->steps[i].score;
    }

    if (total_pos_sum == 0.0) {
        *walkstatpos = 0.0;
        *walkstatneg = 0.0;
        return;
    }

    qsort(ctx->steps, k, sizeof(gene_step_t), compare_gene_steps);

    double running_pos_sum = 0.0;
    double max_peak = 0.0;
    double min_valley = 0.0;
    double neg_step = 1.0 / (double)(n - k);

    for (int i = 0; i < k; i++) {
        int current_rank = ctx->steps[i].rank;
        double neg_penalty = (current_rank - i) * neg_step;

        double normalized_pos = running_pos_sum / total_pos_sum;
        double current_valley = normalized_pos - neg_penalty;
        if (current_valley < min_valley) {
            min_valley = current_valley;
        }

        running_pos_sum += ctx->steps[i].score;

        normalized_pos = running_pos_sum / total_pos_sum;
        double current_peak = normalized_pos - neg_penalty;
        if (current_peak > max_peak) {
            max_peak = current_peak;
        }
    }

    *walkstatpos = max_peak;
    *walkstatneg = min_valley;
}

void
gsva_rnd_walk_nas(int* gsetidx, int k, int* decordstat, double* symrnkstat, int n,
                  double tau, int na_use, int minsize, double* walkstatpos,
                  double* walkstatneg, int* wna) {
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
}

void
gsva_sparse_rnd_walk(int* gsetidx, int k, int* decordstat, double* symrnkstat, int n,
                     double tau, int nzeros, double zerosymrnkstat,
                     double* walkstatpos, double* walkstatneg) {
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
}
