#ifndef CUDA_TYPES_H
#define CUDA_TYPES_H

#include <cuda_runtime.h>
#include <R.h>

#ifdef USE_FP32
    typedef float gsva_float_t;
    #define GSVA_POW(base, exp) powf(base, exp)
    #define GSVA_FABS(val) fabsf(val)
    #define GSVA_EXP(val) expf(val)
    #define GSVA_ROUND(val) roundf(val)
    #define GSVA_NAN ((gsva_float_t)NAN)
#else
    typedef double gsva_float_t;
    #define GSVA_POW(base, exp) pow(base, exp)
    #define GSVA_FABS(val) fabs(val)
    #define GSVA_EXP(val) exp(val)
    #define GSVA_ROUND(val) round(val)
    #define GSVA_NAN NAN
#endif

#define GSVA_BLOCK_C 32
#define GSVA_CUDA_STREAMS 2
#define GSVA_MAX_GPU_GSET_SIZE (8 * 512)

// Number of gene-set size categories.
#define GSVA_NUM_CATS 11

// The single source of truth for gene-set size categories.
#define GSVA_DENSE_WALK(F)                                                     \
  F(32, 1, 32), F(64, 2, 32), F(128, 4, 32), F(160, 5, 32), F(192, 6, 32),    \
  F(256, 4, 64), F(384, 6, 64), F(768, 6, 128), F(1536, 6, 256),              \
  F(2048, 8, 256), F(4096, 8, 512)

#define GSVA_SPARSE_WALK(F)                                                    \
  F(32, 1, 32), F(64, 2, 32), F(128, 4, 32), F(160, 5, 32), F(192, 6, 32),    \
  F(256, 2, 128), F(384, 6, 64), F(768, 6, 128), F(1536, 6, 256),             \
  F(2048, 8, 256), F(4096, 8, 512)

// Accessors to extract single fields from the definitions above.
#define GSVA_CAT_CAP_ENTRY(cap, items, threads) (cap)
#define GSVA_CAT_ITEMS_ENTRY(cap, items, threads) (items)
#define GSVA_CAT_THREADS_ENTRY(cap, items, threads) (threads)

#define GSVA_CUDA_CALL(call) do { \
    cudaError_t err = (call); \
    if (err != cudaSuccess) { \
        error("CUDA error at %s:%d: %s", __FILE__, __LINE__, \
              cudaGetErrorString(err)); \
    } \
} while(0)

#endif
