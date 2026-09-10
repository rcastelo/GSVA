#ifndef RND_WALK_H
#define RND_WALK_H

#ifdef __cplusplus
extern "C" {
#endif

void scatter_walk_single(int* decordstat, double* symrnkstat,
                         int* gset_indices, int gset_size,
                         int N, double tau,
                         double* step_in, double* step_out);

void scan_single(double* step_in, double* step_out, int N);

void compute_walkstats_single(double* step_in, double* step_out,
                              int N, double* walkstats);

void reduce_single(double* walkstats, int N, double* pos, double* neg);

typedef struct {
    int rank;
    double score;
} gene_step_t;

int compare_gene_steps(const void *a, const void *b);

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
