#ifndef CUDA_CBIND_H
#define CUDA_CBIND_H

#include <R.h>
#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif


SEXP
gsva_score_genesets_gpu_R(SEXP ranksR, SEXP genesetsidxR, SEXP intrnksR,
                          SEXP sparseR, SEXP maxdiffR, SEXP absrnkR, SEXP tauR,
                          SEXP minsizeR, SEXP verboseR);

#ifdef __cplusplus
}
#endif

#endif
