#ifndef CUDA_CONFIG_H
#define CUDA_CONFIG_H

#include <math.h>

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

#define GSVA_MAX_GPU_GSET_SIZE (8 * 512)

#endif
