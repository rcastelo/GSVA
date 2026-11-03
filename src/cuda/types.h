#ifndef CUDA_TYPES_H
#define CUDA_TYPES_H

#include <cuda_runtime.h>
#include <R.h>

#ifdef USE_FP32
    typedef float gsva_float_t;
    #define GSVA_POW(base, exp) powf(base, exp)
    #define GSVA_FABS(val) fabsf(val)
    #define GSVA_EXP(val) expf(val)
    #define GSVA_NAN ((gsva_float_t)NAN)
#else
    typedef double gsva_float_t;
    #define GSVA_POW(base, exp) pow(base, exp)
    #define GSVA_FABS(val) fabs(val)
    #define GSVA_EXP(val) exp(val)
    #define GSVA_NAN NAN
#endif

#define GSVA_BLOCK_C 32
#define GSVA_CUDA_STREAMS 2
#define GSVA_MAX_GPU_GSET_SIZE (8 * 512)

#define GSVA_CUDA_CALL(call) do { \
    cudaError_t err = (call); \
    if (err != cudaSuccess) { \
        error("CUDA error at %s:%d: %s", __FILE__, __LINE__, \
              cudaGetErrorString(err)); \
    } \
} while(0)

#endif
