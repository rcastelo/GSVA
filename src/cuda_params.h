#ifndef CUDA_PARAMS_H
#define CUDA_PARAMS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <R.h>
#include <Rinternals.h>

// Device-independent parameters of the CUDA implementation.
SEXP gsva_cuda_max_gset_size_R(void);
SEXP gsva_cuda_walk_safe_R(SEXP GR, SEXP ksR, SEXP tauR);

#ifdef __cplusplus
}
#endif

#endif
