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
 * computational cost from O(p) to O(k). the caller should check beforehand
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

/* random walk on a column of a sparse matrix, where the 'nzeros' genes with a
 * zero expression value share the same rank. genes in the gene set with a
 * nonzero expression value keep their positions in 'decordstat', which lie
 * among the first 'nnz' ones, and those with a zero expression value, all of
 * them with the same weight given by 'zerosymrnkstat', are redistributed
 * evenly across the last 'nzeros' positions of the walk, at the 1-based
 * positions nnz + round(i * step + 1), with step = (nzeros - 1) / (gset_nzeros - 1), for
 * i = 0, ..., gset_nzeros - 1, so that the walk does not depend on the order
 * in which ties among zeros are broken, see .redistribute_geneset_zeros() in
 * R/gsva.R. these positions are distinct, because step >= 1 when
 * gset_nzeros > 1, and are generated in increasing order after the nonzero
 * ones, so only the positions of the nonzero genes need to be sorted. the
 * walk thus has 'k' steps up at distinct positions and p - k steps down,
 * as in the dense random walk, and rnd_walk_extremes() computes its maximum
 * and minimum deviations. sorting only the 'gset_nnz' positions of the
 * nonzero genes in the gene set gives a cost of O(gset_nnz log gset_nnz + k)
 * time instead of O(p), which is at most the O(k log k) of the dense walk.
 *
 * returned value: 1 if the random walk is undefined, 0 otherwise, see
 * gsva_rnd_walk()
 */
int
gsva_sparse_rnd_walk(rndwalk_ctx_t* ctx, int* gsetidx, int k, int* decordstat,
                     double* symrnkstat, int p, double tau, int nzeros,
                     double zerosymrnkstat, double* walkstatpos,
                     double* walkstatneg) {
  int*    ranks = ctx->ranks;
  int*    idx = ctx->idx;
  double* weights = ctx->weights;
  int     nnz = p - nzeros;
  int     gset_nnz = 0;
  int     gset_nzeros;
  double  total_pos_sum = 0.0;

  *walkstatpos = *walkstatneg = NA_REAL;
  if (k == 0 || k == p) /* no in-gene-set or no out-gene-set steps */
    return 1;

  /* genes in the gene set with a nonzero expression value go first */
  for (int i=0; i < k; i++) {
    /* convert 1-based gene indices and positions to 0-based ! */
    int gene_idx = gsetidx[i] - 1;
    int pos = decordstat[gene_idx] - 1;

    if (pos < nnz) { /* assuming nonzero expression values are positive */
      ranks[gset_nnz] = pos;
      idx[gset_nnz] = gset_nnz; /* see gsva_rnd_walk() */
      weights[gset_nnz] = tau == 1 ? symrnkstat[gene_idx] :
                                     pow(symrnkstat[gene_idx], tau);
      total_pos_sum += weights[gset_nnz];
      gset_nnz++;
    }
  }
  gset_nzeros = k - gset_nnz;

  /* R_qsort_int_I() takes 1-based inclusive bounds, so 1, gset_nnz sorts
   * ranks[0..gset_nnz-1] */
  if (gset_nnz > 1)
    R_qsort_int_I(ranks, idx, 1, gset_nnz);

  /* genes in the gene set with a zero expression value go next, at their
   * redistributed positions, which are already in increasing order */
  if (gset_nzeros > 0) {
    double step = 0;
    double zeroweight = tau == 1 ? zerosymrnkstat : pow(zerosymrnkstat, tau);

    if (gset_nzeros > 1)
      step = (double) (nzeros - 1) / (double) (gset_nzeros - 1);

    for (int i=0; i < gset_nzeros; i++) {
      int j = gset_nnz + i;

      /* 1-based position nnz + round(i * step + 1), converted to 0-based */
      ranks[j] = nnz + ((int) round(((double) i) * step + 1)) - 1;
      idx[j] = j;
      weights[j] = zeroweight;
      total_pos_sum += zeroweight;
    }
  }

  if (total_pos_sum == 0.0) /* in-gene-set steps sum up to zero */
    return 1;

  rnd_walk_extremes(ranks, idx, weights, k, p, total_pos_sum, walkstatpos,
                    walkstatneg);

  return 0;
}

/* random walk on a column with missing values. 'p' is the number of
 * nonmissing genes in the column, i.e., the length of the random walk, since
 * nonmissing genes have decreasing order statistics between 1 and 'p', and
 * missing genes are left out of the walk. with na_use < 3 (everything or
 * all.obs), a gene set with missing values gives NA, and with na_use == 3
 * (na.rm), its missing genes are skipped and it gives NA, and sets '*wna', if
 * it has less than 'minsize' nonmissing genes. otherwise, the 'k_notna'
 * nonmissing genes of the gene set are sorted by position, as in
 * gsva_rnd_walk(), and rnd_walk_extremes() computes the maximum and minimum
 * deviations of the walk in O(k log k) time instead of O(p).
 *
 * returned value: 1 if the random walk is undefined, 0 otherwise, see
 * gsva_rnd_walk(); a gene set with missing values under na_use < 3 or with
 * less than 'minsize' nonmissing genes does not lead to an undefined walk
 */
int
gsva_rnd_walk_nas(rndwalk_ctx_t* ctx, int* gsetidx, int k, int* decordstat,
                  double* symrnkstat, int p, double tau, int na_use,
                  int minsize, double* walkstatpos, double* walkstatneg,
                  int* wna) {
  int*    ranks = ctx->ranks;
  int*    idx = ctx->idx;
  double* weights = ctx->weights;
  int     k_notna = 0;
  double  total_pos_sum = 0.0;

  *walkstatpos = *walkstatneg = NA_REAL;

  for (int i=0; i < k; i++) {
    /* convert 1-based gene indices and positions to 0-based ! */
    int gene_idx = gsetidx[i] - 1;

    if (decordstat[gene_idx] != NA_INTEGER) { /* na.rm skips NAs */
      ranks[k_notna] = decordstat[gene_idx] - 1;
      idx[k_notna] = k_notna; /* see gsva_rnd_walk() */
      weights[k_notna] = tau == 1 ? symrnkstat[gene_idx] :
                                    pow(symrnkstat[gene_idx], tau);
      total_pos_sum += weights[k_notna];
      k_notna++;
    } else if (na_use < 3) /* everything or all.obs */
      return 0;
  }

  if (k_notna < minsize) { /* na.rm */
    *wna = 1;
    return 0;
  }

  /* no in-gene-set or no out-gene-set steps, or in-gene-set steps sum up
   * to zero */
  if (k_notna == 0 || k_notna == p || total_pos_sum == 0.0)
    return 1;

  /* R_qsort_int_I() takes 1-based inclusive bounds, so 1, k_notna sorts
   * ranks[0..k_notna-1] */
  R_qsort_int_I(ranks, idx, 1, k_notna);

  rnd_walk_extremes(ranks, idx, weights, k_notna, p, total_pos_sum,
                    walkstatpos, walkstatneg);

  return 0;
}
