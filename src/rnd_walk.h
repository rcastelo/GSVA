#ifndef RND_WALK_H
#define RND_WALK_H

#include <R.h>
#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int rank;
    double score;
} gene_step_t;

typedef struct {
  int         maxk;   /* size of the largest gene set */
  gene_step_t *steps; /* preallocated, capacity maxk */
} rndwalk_ctx_t;

rndwalk_ctx_t* rndwalk_ctx_create(SEXP genesetsidxR);

void gsva_rnd_walk(rndwalk_ctx_t* ctx, int* gsetidx, int k,
                   int* decordstat, double* symrnkstat, int n,
                   double tau, double* walkstatpos, double* walkstatneg);

void gsva_rnd_walk_nas(int* gsetidx, int k, int* decordstat, double* symrnkstat, int n,
                       double tau, int na_use, int minsize, double* walkstatpos,
                       double* walkstatneg, int* wna);
void gsva_sparse_rnd_walk(int* gsetidx, int k, int* decordstat, double* symrnkstat, int n,
                          double tau, int nzeros, double zerosymrnkstat,
                          double* walkstatpos, double* walkstatneg);

#ifdef __cplusplus
}
#endif

#endif
