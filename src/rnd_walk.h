#ifndef RND_WALK_H
#define RND_WALK_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int rank;
    double score;
} gene_step_t;


void gsva_rnd_walk(int* gsetidx, int k, int* decordstat, double* symrnkstat, int n,
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
