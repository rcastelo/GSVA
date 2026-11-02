#ifndef RND_WALK_H
#define RND_WALK_H

#include <R.h>
#include <Rinternals.h>

/* context for the dense and sparse random walks, created once per call to
 * gsva_score_genesets_R(), with a buffer sized for the largest gene set
 * that is reused for every gene set and column */
typedef struct {
  int     maxk;    /* size of the largest gene set */
  int*    ranks;   /* 0-based positions of the genes in the gene set */
  int*    idx;     /* permutation applied to 'ranks' by the sort */
  double* weights; /* weights of the genes in the gene set */
} rndwalk_ctx_t;

rndwalk_ctx_t*
rndwalk_ctx_create(SEXP genesetsidxR);

int
gsva_rnd_walk(rndwalk_ctx_t* ctx, int* gsetidx, int k, int* decordstat,
              double* symrnkstat, int p, double tau, double* walkstatpos,
              double* walkstatneg);

int
gsva_sparse_rnd_walk(rndwalk_ctx_t* ctx, int* gsetidx, int k, int* decordstat,
                     double* symrnkstat, int p, double tau, int nzeros,
                     double zerosymrnkstat, double* walkstatpos,
                     double* walkstatneg);

int
gsva_rnd_walk_nas(int* gsetidx, int k, int* decordstat, double* symrnkstat,
                  int n, double tau, int na_use, int minsize,
                  double* walkstatpos, double* walkstatneg, int* wna);

#endif
