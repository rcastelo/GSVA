#include "cuda/cbind.h"
#include "cuda/kernels.cuh"
#include "cuda/memory.h"
#include "cuda/types.h"
#include "ranks.h"
#include "rnd_walk.h"
#include <R.h>
#include <Rdefines.h>
#include <math.h>
#include <cuda_runtime.h>
#include <time.h>

static inline int min2(int a, int b) { return a < b ? a : b; }
static inline int max2(int a, int b) { return a > b ? a : b; }

/* Pack sparse DGCMatrix columns into COO-style (gene_offset, value) arrays.
 * Returns total_nnz across the block. */
static int pack_ranks_dgc(ranks_ctx_t *ctx, gsva_host_t *host,
                          int stream_idx, int block_c, int block_size) {
  const int *dgc_i = ctx->u.dgc.i;
  const int *dgc_p = ctx->u.dgc.p;
  const double *dgc_x = ctx->u.dgc.x;
  int *out_offt = host->sparse_offt[stream_idx];
  int *out_vals = host->sparse_vals[stream_idx];
  int *out_cell_offt = host->cell_offt[stream_idx];
  int *out_nnzpercell = host->nnzpercell[stream_idx];

  int offset = 0;
  for (int c = 0; c < block_size; c++) {
    int col = block_c + c;
    out_cell_offt[c] = offset;
    out_nnzpercell[c] = dgc_p[col + 1] - dgc_p[col];
    for (int idx = dgc_p[col]; idx < dgc_p[col + 1]; idx++) {
      out_offt[offset] = dgc_i[idx];
      out_vals[offset] = (int)dgc_x[idx];
      offset++;
    }
  }
  out_cell_offt[block_size] = offset;
  return offset;
}

/* Pack sparse SVT columns into COO-style (gene_offset, value) arrays.
 * Returns total_nnz across the block. Whether vals come from REAL (double)
 * or INTEGER is inferred from ctx->type. */
static int pack_ranks_svt(ranks_ctx_t *ctx, gsva_host_t *host,
                          int stream_idx, int block_c, int block_size) {
  SEXP svt = ctx->u.svt;
  Rboolean is_dbl = (Rboolean)(ctx->type == RANKSTYPE_SVT_DBL);
  int *out_offt = host->sparse_offt[stream_idx];
  int *out_vals = host->sparse_vals[stream_idx];
  int *out_cell_offt = host->cell_offt[stream_idx];
  int *out_nnzpercell = host->nnzpercell[stream_idx];

  int offset = 0;
  for (int c = 0; c < block_size; c++) {
    int col = block_c + c;
    SEXP leaf = VECTOR_ELT(svt, col);
    out_cell_offt[c] = offset;
    if (leaf == R_NilValue) {
      out_nnzpercell[c] = 0;
    } else {
      SEXP valsR = VECTOR_ELT(leaf, 0);
      SEXP offsetsR = VECTOR_ELT(leaf, 1);
      int nvals = length(valsR);
      int noffsets = length(offsetsR);
      int *offsets = INTEGER(offsetsR);
      if (nvals > 0) {
        if (is_dbl) {
          double *vals_dbl = REAL(valsR);
          for (int k = 0; k < nvals; k++) {
            out_offt[offset] = offsets[k];
            out_vals[offset] = (int)vals_dbl[k];
            offset++;
          }
        } else {
          int *vals = INTEGER(valsR);
          for (int k = 0; k < nvals; k++) {
            out_offt[offset] = offsets[k];
            out_vals[offset] = vals[k];
            offset++;
          }
        }
      } else {
        /* lacunar: all values = 1 */
        for (int k = 0; k < noffsets; k++) {
          out_offt[offset] = offsets[k];
          out_vals[offset] = 1;
          offset++;
        }
      }
      out_nnzpercell[c] = noffsets;
    }
  }
  out_cell_offt[block_size] = offset;
  return offset;
}

static void ranks2stats_gpu(ranks_ctx_t *ctx, int G, int block_c,
                            int block_size, int sparse_mode,
                            gsva_device_t *device, gsva_host_t *host,
                            int stream_idx, SEXP ranksR) {
  cudaStream_t stream = device->stream[stream_idx];
  r2s_fetch_type_t fetch_type;

  int *h_offt = host->sparse_offt[stream_idx];
  int *h_vals = host->sparse_vals[stream_idx];
  int *h_cofft = host->cell_offt[stream_idx];
  int *h_nnzc = host->nnzpercell[stream_idx];

  size_t sz_cofft = (size_t)(block_size + 1) * sizeof(int);
  size_t sz_nnzc = (size_t)block_size * sizeof(int);

  if (ctx->type == RANKSTYPE_MATRIX_INT) {
    /* Dense integer matrix: transfer from column-major R to row-major device */
    fetch_type = R2S_DENSE;
    const int *ranks_data = INTEGER(ranksR);
    for (int c = 0; c < block_size; c++) {
      GSVA_CUDA_CALL(cudaMemcpyAsync(device->dense_ranks[stream_idx] + c * G,
                                     ranks_data + (block_c + c) * G,
                                     (size_t)G * sizeof(int),
                                     cudaMemcpyHostToDevice, stream));
    }
  } else if (ctx->type == RANKSTYPE_MATRIX_DBL) {
    /* Dense double matrix: cast double to int for all cells */
    fetch_type = R2S_DENSE;
    const double *ranks_data_dbl = REAL(ranksR);
    int *h_cast = host->dense_cast[stream_idx];  // [G * block_size]
    for (int c = 0; c < block_size; c++) {
      const double *src = ranks_data_dbl + (block_c + c) * G;
      for (int g = 0; g < G; g++) {
        h_cast[c * G + g] = (int)src[g];
      }
    }
    GSVA_CUDA_CALL(cudaMemcpyAsync(device->dense_ranks[stream_idx], h_cast,
                                   (size_t)block_size * G * sizeof(int),
                                   cudaMemcpyHostToDevice, stream));
  } else if (ctx->type == RANKSTYPE_DGC) {
    fetch_type = R2S_SPARSE;
    int total_nnz =
        pack_ranks_dgc(ctx, host, stream_idx, block_c, block_size);
    GSVA_CUDA_CALL(cudaMemcpyAsync(device->sparse_offt[stream_idx], h_offt,
                                   (size_t)total_nnz * sizeof(int),
                                   cudaMemcpyHostToDevice, stream));
    GSVA_CUDA_CALL(cudaMemcpyAsync(device->sparse_vals[stream_idx], h_vals,
                                   (size_t)total_nnz * sizeof(int),
                                   cudaMemcpyHostToDevice, stream));
    GSVA_CUDA_CALL(cudaMemcpyAsync(device->cell_offt[stream_idx], h_cofft,
                                   sz_cofft, cudaMemcpyHostToDevice, stream));
    GSVA_CUDA_CALL(cudaMemcpyAsync(device->nnzpercell[stream_idx], h_nnzc,
                                   sz_nnzc, cudaMemcpyHostToDevice, stream));
  } else if (ctx->type == RANKSTYPE_SVT_INT || ctx->type == RANKSTYPE_SVT_DBL) {
    fetch_type = R2S_SPARSE;
    int total_nnz = pack_ranks_svt(ctx, host, stream_idx, block_c, block_size);
    GSVA_CUDA_CALL(cudaMemcpyAsync(device->sparse_offt[stream_idx], h_offt,
                                   (size_t)total_nnz * sizeof(int),
                                   cudaMemcpyHostToDevice, stream));
    GSVA_CUDA_CALL(cudaMemcpyAsync(device->sparse_vals[stream_idx], h_vals,
                                   (size_t)total_nnz * sizeof(int),
                                   cudaMemcpyHostToDevice, stream));
    GSVA_CUDA_CALL(cudaMemcpyAsync(device->cell_offt[stream_idx], h_cofft,
                                   sz_cofft, cudaMemcpyHostToDevice, stream));
    GSVA_CUDA_CALL(cudaMemcpyAsync(device->nnzpercell[stream_idx], h_nnzc,
                                   sz_nnzc, cudaMemcpyHostToDevice, stream));
  } else {
    error("Unknown ranks_ctx type: %d", ctx->type);
  }

  int threads = (fetch_type == R2S_SPARSE) ? GSVA_R2S_THREADS : min2(G, GSVA_R2S_THREADS);
  ranks2stats_gpu<<<block_size, threads, 0, stream>>>(
      fetch_type, device->dense_ranks[stream_idx],
      device->sparse_offt[stream_idx], device->sparse_vals[stream_idx],
      device->cell_offt[stream_idx], device->nnzpercell[stream_idx],
      device->r_scratch[stream_idx], device->decordstat[stream_idx],
      device->symrnkstat[stream_idx], G, block_size, (int)(sparse_mode != 0));

  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    error("ranks2stats_gpu launch error at %s:%d: %s", __FILE__, __LINE__,
          cudaGetErrorString(err));
  }
}

static void gsva_rnd_walk_gpu_dispatch(gsva_float_t *h_es, int S, int G, int C,
                                       gsva_float_t tau, int score_type,
                                       gsva_device_t *device, int stream_idx) {
  cudaStream_t stream = device->stream[stream_idx];
  size_t es_size = (size_t)S * C * sizeof(gsva_float_t);
  int max_k = device->max_k;

  dim3 numBlocks(C, S);
  cudaError_t err;

  if (max_k <= 32) {
    gsea_walk_kernel<1, 32><<<numBlocks, 32, 0, stream>>>(
        device->gsetofft, device->gsetidxs, device->decordstat[stream_idx],
        device->symrnkstat[stream_idx], device->es[stream_idx], S, G, C, tau,
        score_type);
  } else if (max_k <= 128) {
    gsea_walk_kernel<4, 32><<<numBlocks, 32, 0, stream>>>(
        device->gsetofft, device->gsetidxs, device->decordstat[stream_idx],
        device->symrnkstat[stream_idx], device->es[stream_idx], S, G, C, tau,
        score_type);
  } else if (max_k <= 2048) {
    gsea_walk_kernel<4, 512><<<numBlocks, 512, 0, stream>>>(
        device->gsetofft, device->gsetidxs, device->decordstat[stream_idx],
        device->symrnkstat[stream_idx], device->es[stream_idx], S, G, C, tau,
        score_type);
  } else {
    gsea_walk_kernel<8, 512><<<numBlocks, 512, 0, stream>>>(
        device->gsetofft, device->gsetidxs, device->decordstat[stream_idx],
        device->symrnkstat[stream_idx], device->es[stream_idx], S, G, C, tau,
        score_type);
  }

  err = cudaGetLastError();
  if (err != cudaSuccess) {
    error("CUDA kernel launch error: %s", cudaGetErrorString(err));
  }

  GSVA_CUDA_CALL(cudaMemcpyAsync(h_es, device->es[stream_idx], es_size,
                                 cudaMemcpyDeviceToHost, stream));
}

typedef struct {
  SEXP ranksR;
  SEXP genesetsidxR;
  Rboolean intrnks;
  Rboolean sparse;
  gsva_float_t tau;
  int score_type;
  gsva_device_t device;
  gsva_host_t host;
} gsva_job_t;

static SEXP gsva_score_body(void *data) {
  gsva_job_t *job = (gsva_job_t *)data;

  ranks_ctx_t *ctx = ranks_ctx_create(job->ranksR, job->intrnks, job->sparse);
  int G = ctx->p;
  int C = ctx->n;

  gsva_device_init(&job->device, job->genesetsidxR, G, job->sparse);
  int S = job->device.S;
  gsva_host_init(&job->host, G, S, job->sparse);

  SEXP esR;
  PROTECT(esR = allocMatrix(REALSXP, S, C));
  double *es = REAL(esR);

  int total_blocks = (C + GSVA_BLOCK_C - 1) / GSVA_BLOCK_C;
  for (int block_id = 0; block_id < total_blocks; block_id++) {
    int stream_idx = block_id % GSVA_CUDA_STREAMS;
    int block_c = block_id * GSVA_BLOCK_C;
    int block_c_end = min2(block_c + GSVA_BLOCK_C, C);
    int block_size = block_c_end - block_c;

    GSVA_CUDA_CALL(cudaStreamSynchronize(job->device.stream[stream_idx]));

    if (block_id >= GSVA_CUDA_STREAMS) {
      int prev_block = block_id - GSVA_CUDA_STREAMS;
      int prev_c = prev_block * GSVA_BLOCK_C;
      int prev_end = min2(prev_c + GSVA_BLOCK_C, C);
      int prev_size = prev_end - prev_c;

      gsva_float_t *prev_es = job->host.es[stream_idx];
      for (int c = 0; c < prev_size; c++) {
        for (int s = 0; s < S; s++) {
          double v = (double)prev_es[c * S + s];
          if (isnan(v)) {
            // Convert to a R-level NA
            es[(prev_c + c) * S + s] = NA_REAL;
          } else {
            es[(prev_c + c) * S + s] = v;
          }
        }
      }
    }

    ranks2stats_gpu(ctx, G, block_c, block_size, (int)job->sparse,
                    &job->device, &job->host, stream_idx, job->ranksR);
    gsva_rnd_walk_gpu_dispatch(job->host.es[stream_idx], S, G, block_size,
                               job->tau, job->score_type, &job->device,
                               stream_idx);
  }

  cudaDeviceSynchronize();
  for (int block_id = max2(0, total_blocks - GSVA_CUDA_STREAMS);
       block_id < total_blocks; block_id++) {
    int stream_idx = block_id % GSVA_CUDA_STREAMS;
    int block_c = block_id * GSVA_BLOCK_C;
    int block_c_end = min2(block_c + GSVA_BLOCK_C, C);
    int block_size = block_c_end - block_c;

    gsva_float_t *h_es = job->host.es[stream_idx];
    for (int c = 0; c < block_size; c++) {
      for (int s = 0; s < S; s++) {
        double v = (double)h_es[c * S + s];
        if (isnan(v)) {
          // Convert to a R-level NA
          es[(block_c + c) * S + s] = NA_REAL;
        } else {
          es[(block_c + c) * S + s] = v;
        }
      }
    }
  }

  UNPROTECT(1);
  return esR;
}

static void gsva_score_cleanup(void *data) {
  gsva_job_t *job = (gsva_job_t *)data;
  gsva_host_deinit(&job->host);
  gsva_device_deinit(&job->device);
}

SEXP gsva_score_genesets_gpu_R(SEXP ranksR, SEXP genesetsidxR, SEXP intrnksR,
                               SEXP sparseR, SEXP maxdiffR, SEXP absrnkR,
                               SEXP tauR, SEXP minsizeR, SEXP verboseR) {
  gsva_job_t job = {
      .ranksR = ranksR,
      .genesetsidxR = genesetsidxR,
      .intrnks = (Rboolean)asLogical(intrnksR),
      .sparse = (Rboolean)asLogical(sparseR),
      .tau = (gsva_float_t)REAL(tauR)[0],
      .score_type = (Rboolean)asLogical(maxdiffR)
                        ? ((Rboolean)asLogical(absrnkR) ? 1 : 0)
                        : 2
  };

  return R_ExecWithCleanup(gsva_score_body, &job, gsva_score_cleanup, &job);
}
