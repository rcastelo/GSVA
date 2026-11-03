#ifndef CUDA_MEMORY_H
#define CUDA_MEMORY_H

#include "cuda/types.h"
#include <R.h>
#include <Rdefines.h>
#include <cuda.h>
#include <cuda_runtime.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  int G, S, max_k;
  int *gsetofft, *gsetidxs;

  cudaStream_t stream[GSVA_CUDA_STREAMS];

  int *r_scratch[GSVA_CUDA_STREAMS];   // [G*B_C] scratch for rank shift

  // Buffers for sparse ranks
  int *sparse_offt[GSVA_CUDA_STREAMS]; // [G*B_C] gene offsets
  int *sparse_vals[GSVA_CUDA_STREAMS]; // [G*B_C] rank values
  int *cell_offt[GSVA_CUDA_STREAMS];   // [B_C+1] prefix sum of nnz
  int *nnzpercell[GSVA_CUDA_STREAMS];  // [B_C] nnz per cell

  // Buffer for dense ranks
  int *dense_ranks[GSVA_CUDA_STREAMS]; // [G*B_C] dense ranks

  // Buffers for random walk
  int *decordstat[GSVA_CUDA_STREAMS];  // [G*B_C] decreasing order statistics
  gsva_float_t
      *symrnkstat[GSVA_CUDA_STREAMS];  // [G*B_C] symmetric rank statistics
  gsva_float_t *es[GSVA_CUDA_STREAMS]; // [S*B_C] enrichment scores
} gsva_device_t;

void gsva_device_init(gsva_device_t *d, SEXP genesetsidxR, int G,
                      Rboolean sparse);
void gsva_device_deinit(gsva_device_t *d);

typedef struct {
  int *sparse_offt[GSVA_CUDA_STREAMS]; // [G*B_C] gene offsets
  int *sparse_vals[GSVA_CUDA_STREAMS]; // [G*B_C] rank values
  int *cell_offt[GSVA_CUDA_STREAMS];   // [B_C+1] prefix sum
  int *nnzpercell[GSVA_CUDA_STREAMS];  // [B_C] nnz per cell

  int *dense_cast[GSVA_CUDA_STREAMS];  // [G*B_C] workspace for double-to-int cast

  gsva_float_t *es[GSVA_CUDA_STREAMS]; // [S*B_C] enrichment scores
} gsva_host_t;

void gsva_host_init(gsva_host_t *h, int G, int S, Rboolean sparse);
void gsva_host_deinit(gsva_host_t *h);

#ifdef __cplusplus
}
#endif

#endif
