#ifndef KERNELS_CUH
#define KERNELS_CUH

#include <math.h>
#include <limits.h>
#include <cub/cub.cuh>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include "cuda/types.h"

#define GSVA_R2S_THREADS 512

typedef enum {
    R2S_DENSE  = 1,  // Dense rank matrix (int or double)
    R2S_SPARSE = 2   // Sparse layout: DGC or SVT_INT/DGC
} r2s_fetch_type_t;

__global__ void ranks2stats_gpu(
    r2s_fetch_type_t fetch_type,
    const int* __restrict__ dense_ranks,
    const int* __restrict__ sparse_offt,
    const int* __restrict__ sparse_vals,
    const int* __restrict__ cell_offt,
    const int* __restrict__ nnzpercell,
    int* __restrict__ d_r,
    int* __restrict__ decordstat,
    gsva_float_t* __restrict__ symrnkstat,
    int G,
    int block_size,
    int sparse_mode
) {
    int cell = blockIdx.x;
    int tid  = threadIdx.x;
    int threads = blockDim.x;
    if (cell >= block_size) return;

    int nnz, nzs;
    if (fetch_type == R2S_DENSE) {
        nnz = G;
        nzs = 0;
        const int* dr = dense_ranks + cell * G;
        for (int g = tid; g < G; g += threads) {
            d_r[cell * G + g] = dr[g];
        }
    } else {
        int base = cell_offt[cell];
        nnz = nnzpercell[cell];
        nzs = G - nnz;

        for (int g = tid; g < G; g += threads) {
            d_r[cell * G + g] = 0;
        }
        __syncthreads();

        for (int i = tid; i < nnz; i += threads) {
            int pos    = base + i;
            int gene   = sparse_offt[pos];
            int val    = sparse_vals[pos];
            d_r[cell * G + gene] = val;
        }
    }
    __syncthreads();

    if (fetch_type == R2S_SPARSE) {
        struct ZeroCountPrefixOp {
            int running_total;
            __device__ ZeroCountPrefixOp(int r) : running_total(r) {}
            __device__ int operator()(int chunk_zeros) {
                int old = running_total;
                running_total += chunk_zeros;
                return old;
            }
        };

        ZeroCountPrefixOp prefix_op(0);
        int num_chunks = (G + threads - 1) / threads;

        using BlockScan = cub::BlockScan<int, GSVA_R2S_THREADS>;
        __shared__ typename BlockScan::TempStorage scan_storage;

        for (int chunk = 0; chunk < num_chunks; chunk++) {
            int idx = chunk * threads + tid;

            int cur_val = (idx < G) ? d_r[cell * G + idx] : 0;
            int is_zero = (cur_val == 0) ? 1 : 0;

            int scan_result = 0;
            BlockScan(scan_storage).ExclusiveSum(is_zero, scan_result, prefix_op);

            if (idx < G) {
                if (is_zero) {
                    d_r[cell * G + idx] = scan_result + 1;
                } else {
                    d_r[cell * G + idx] = cur_val + nzs;
                }
            }
            __syncthreads();
        }
    }

    for (int g = tid; g < G; g += threads) {
        int rshifted = d_r[cell * G + g];
        decordstat[cell * G + g] = G - rshifted + 1;

        if (sparse_mode && nzs > 0 && fetch_type == R2S_SPARSE) {
            if (rshifted > nzs) {
                int raw_r = rshifted - nzs;
                symrnkstat[cell * G + g] = GSVA_FABS((gsva_float_t)(nnz + 1) / 2.0f - (gsva_float_t)(raw_r + 1));
            } else {
                /* Was zero */
                symrnkstat[cell * G + g] = GSVA_FABS((gsva_float_t)(nnz + 1) / 2.0f - 1.0f);
            }
        } else {
            symrnkstat[cell * G + g] = GSVA_FABS((gsva_float_t)G / 2.0f - (gsva_float_t)rshifted);
        }
    }
}

struct PeakValley {
    gsva_float_t peak;
    gsva_float_t valley;
};

struct GSVAMaxOp {
    __device__ gsva_float_t operator()(gsva_float_t a, gsva_float_t b) const {
        return a > b ? a : b;
    }
};

struct PVReduceOp {
    __device__ PeakValley operator()(PeakValley a, PeakValley b) const {
        return {
            a.peak > b.peak ? a.peak : b.peak,
            a.valley < b.valley ? a.valley : b.valley
        };
    }
};

template<int GSVA_WALK_ITEMS, int GSVA_WALK_THREADS>
__global__ void gsea_dense_walk_kernel(
    const int* __restrict__ gsetofft,
    const int* __restrict__ gsetidxs,
    const int* __restrict__ decordstat_block,
    const gsva_float_t* __restrict__ symrnkstat_block,
    gsva_float_t* __restrict__ out_es,
    const int* __restrict__ cat_gset,
    int S,
    int G,
    int C,
    gsva_float_t tau,
    int score_type
) {
    typedef cub::BlockRadixSort<int, GSVA_WALK_THREADS, GSVA_WALK_ITEMS, gsva_float_t> BlockRadixSort;
    typedef cub::BlockScan<double, GSVA_WALK_THREADS> BlockScanDouble;
    typedef cub::BlockReduce<PeakValley, GSVA_WALK_THREADS> BlockReducePV;

    union TempStorage {
        typename BlockRadixSort::TempStorage sort_storage;
        typename BlockScanDouble::TempStorage scan_storage;
        typename BlockReducePV::TempStorage reduce_storage;
    };

    __shared__ TempStorage tmp;
    __shared__ double s_thread_total;

    int c = blockIdx.x;
    int gset_local = blockIdx.y;
    int tid = threadIdx.x;

    int gset = cat_gset[gset_local];
    if (c >= C || gset_local >= S) return;

    int offset = gsetofft[gset];
    int k = gsetofft[gset + 1] - offset;

    int keys[GSVA_WALK_ITEMS];
    gsva_float_t vals[GSVA_WALK_ITEMS];

    int thread_start = tid * GSVA_WALK_ITEMS;
    for (int i = 0; i < GSVA_WALK_ITEMS; i++) {
        int idx = thread_start + i;
        if (idx < k) {
            int gene_idx = gsetidxs[offset + idx];
            int block_idx = c * G + (gene_idx - 1);
            keys[i] = decordstat_block[block_idx] - 1;
            vals[i] = GSVA_POW(symrnkstat_block[block_idx], tau);
        } else {
            keys[i] = INT_MAX;
            vals[i] = 0.0f;
        }
    }

    BlockRadixSort(tmp.sort_storage).Sort(keys, vals);
    __syncthreads();

    double stat_sums[GSVA_WALK_ITEMS];
    for (int i = 0; i < GSVA_WALK_ITEMS; i++) {
        if (keys[i] != INT_MAX) {
            stat_sums[i] = vals[i];
        } else {
            stat_sums[i] = 0.0;
        }
    }

    BlockScanDouble(tmp.scan_storage).InclusiveSum(stat_sums, stat_sums);
    if (tid == GSVA_WALK_THREADS - 1) s_thread_total = stat_sums[GSVA_WALK_ITEMS - 1];
    __syncthreads();

    double stat_total = s_thread_total;
    if (!(stat_total > 0.0 && isfinite(stat_total)) || k == G) {
        if (tid == 0) out_es[c * S + gset] = GSVA_NAN;
        return;
    }

    double thread_max_peak = -1e30;
    double thread_min_valley = 1e30;
    double neg_step = 1.0f / (gsva_float_t)(G - k);
    double inv_stat_total = 1.0 / stat_total;

    for (int i = 0; i < GSVA_WALK_ITEMS; i++) {
        int idx = thread_start + i;
        if (keys[i] == INT_MAX) continue;

        double stat_running = stat_sums[i] - vals[i];
        double stat_norm = stat_running * inv_stat_total;
        double neg_penalty = (keys[i] - idx) * neg_step;
        double current_valley = stat_norm - neg_penalty;
        if (current_valley < thread_min_valley) thread_min_valley = current_valley;

        stat_running += vals[i];
        stat_norm = stat_running * inv_stat_total;
        double current_peak = stat_norm - neg_penalty;
        if (current_peak > thread_max_peak) thread_max_peak = current_peak;
    }

    PeakValley pv = {(gsva_float_t)thread_max_peak, (gsva_float_t)thread_min_valley};
    PeakValley result = BlockReducePV(tmp.reduce_storage).Reduce(pv, PVReduceOp());

    if (tid == 0) {
        gsva_float_t final_max_peak = result.peak;
        gsva_float_t final_min_valley = result.valley;
        gsva_float_t res = 0.0f;
        switch (score_type) {
            case 0: res = final_max_peak + final_min_valley; break;
            case 1: res = final_max_peak - final_min_valley; break;
            case 2: res = (final_max_peak > GSVA_FABS(final_min_valley)) ? final_max_peak : final_min_valley; break;
            default: res = -1.0f; break;
        }
        out_es[c * S + gset] = res;
    }
}

template<int GSVA_WALK_ITEMS, int GSVA_WALK_THREADS>
__global__ void gsea_sparse_walk_kernel(
    const int* __restrict__ gsetofft,
    const int* __restrict__ gsetidxs,
    const int* __restrict__ decordstat_block,
    const gsva_float_t* __restrict__ symrnkstat_block,
    gsva_float_t* __restrict__ out_es,
    const int* __restrict__ cat_gset,
    int S, int G, int C,
    gsva_float_t tau, int score_type,
    const int* __restrict__ nnzpercell
) {
    typedef cub::BlockRadixSort<int, GSVA_WALK_THREADS, GSVA_WALK_ITEMS, gsva_float_t> BlockRadixSort;
    typedef cub::BlockScan<double, GSVA_WALK_THREADS> BlockScanDouble;
    typedef cub::BlockScan<int, GSVA_WALK_THREADS> BlockScanInt;
    typedef cub::BlockReduce<PeakValley, GSVA_WALK_THREADS> BlockReducePV;
    typedef cub::BlockReduce<int, GSVA_WALK_THREADS> BlockReduceInt;
    typedef cub::BlockReduce<gsva_float_t, GSVA_WALK_THREADS> BlockReduceFloat;

    union TempStorage {
        typename BlockRadixSort::TempStorage sort_storage;
        typename BlockScanDouble::TempStorage scan_float;
        typename BlockScanInt::TempStorage scan_int;
        typename BlockReduceInt::TempStorage reduce_int;
        typename BlockReduceFloat::TempStorage reduce_float;
        typename BlockReducePV::TempStorage reduce;
    };

    __shared__ TempStorage tmp;
    __shared__ double s_nonzero_total;
    __shared__ int s_nonzero_events;
    __shared__ int s_gset_nonzeros;
    __shared__ gsva_float_t s_zero_weight;

    int c = blockIdx.x;
    int gset_local = blockIdx.y;
    int tid = threadIdx.x;
    if (c >= C || gset_local >= S) return;

    int gset = cat_gset[gset_local];
    int offset = gsetofft[gset];
    int k = gsetofft[gset + 1] - offset;

    int nnz = nnzpercell[c];
    int nzeros = G - nnz;
    int base = c * G;

    int keys[GSVA_WALK_ITEMS];
    gsva_float_t vals[GSVA_WALK_ITEMS];
    int zcount = 0;
    gsva_float_t zero_weight = (gsva_float_t)0.0;

    int thread_start = tid * GSVA_WALK_ITEMS;
    for (int i = 0; i < GSVA_WALK_ITEMS; i++) {
        int idx = thread_start + i;
        if (idx < k) {
            int gene = gsetidxs[offset + idx] - 1;
            int cell = base + gene;
            int desc_rank = decordstat_block[cell];  // 1-based
            if (desc_rank <= nnz) {
                keys[i] = desc_rank - 1;  // 0-based
                vals[i] = GSVA_POW(symrnkstat_block[cell], tau);
            } else {
                // in-set, zero-expressed gene: kept out of the sort array
                zcount++;
                if (zcount == 1)
                    zero_weight = GSVA_POW(symrnkstat_block[cell], tau);
                
                keys[i] = INT_MAX;
                vals[i] = (gsva_float_t)0.0;
            }
        } else {
            keys[i] = INT_MAX;
            vals[i] = 0.0f;
        }
    }

    int gset_nonzeros = BlockReduceInt(tmp.reduce_int).Sum(zcount);
    if (tid == 0) s_gset_nonzeros = gset_nonzeros;
    __syncthreads();
    gset_nonzeros = s_gset_nonzeros;

    zero_weight = BlockReduceFloat(tmp.reduce_float).Reduce(zero_weight, GSVAMaxOp());
    if (tid == 0) s_zero_weight = zero_weight;
    __syncthreads();
    zero_weight = s_zero_weight;

    BlockRadixSort(tmp.sort_storage).Sort(keys, vals);
    __syncthreads();

    int keep[GSVA_WALK_ITEMS];
    double stat_sums[GSVA_WALK_ITEMS];
    for (int i = 0; i < GSVA_WALK_ITEMS; i++) {
        keep[i] = (keys[i] != INT_MAX) ? 1 : 0;
        stat_sums[i] = keep[i] ? vals[i] : 0.0;
    }

    int event_idxs[GSVA_WALK_ITEMS];
    BlockScanInt(tmp.scan_int).InclusiveSum(keep, event_idxs);
    if (tid == GSVA_WALK_THREADS - 1) s_nonzero_events = event_idxs[GSVA_WALK_ITEMS - 1];
    __syncthreads();
    int nonzero_events = s_nonzero_events;

    BlockScanDouble(tmp.scan_float).InclusiveSum(stat_sums, stat_sums);
    if (tid == GSVA_WALK_THREADS - 1) s_nonzero_total = stat_sums[GSVA_WALK_ITEMS - 1];
    __syncthreads();

    double nonzero_total = s_nonzero_total;
    double stat_total = nonzero_total + (double)gset_nonzeros * zero_weight;
    if (!(stat_total > 0.0 && isfinite(stat_total)) || k == G) {
        if (tid == 0) out_es[c * S + gset] = GSVA_NAN;
        return;
    }

    double thread_max_peak = -1e30;
    double thread_min_valley = 1e30;
    gsva_float_t neg_step = 1.0f / (gsva_float_t)(G - (nonzero_events + gset_nonzeros));
    double inv_stat_total = 1.0 / stat_total;
    double step = 0.0;
    if (gset_nonzeros > 1)
        step = ((double) nzeros - 1.0) / ((double) gset_nonzeros - 1.0);

    // in-set, non-zero expressed genes
    for (int i = 0; i < GSVA_WALK_ITEMS; i++) {
        if (!keep[i]) continue;
        int idx = event_idxs[i] - 1;

        double stat_running = stat_sums[i] - vals[i];
        double stat_norm = stat_running * inv_stat_total;
        gsva_float_t neg_penalty = (keys[i] - idx) * neg_step;
        double current_valley = stat_norm - neg_penalty;
        if (current_valley < thread_min_valley) thread_min_valley = current_valley;

        stat_running += vals[i];
        stat_norm = stat_running * inv_stat_total;
        double current_peak = stat_norm - neg_penalty;
        if (current_peak > thread_max_peak) thread_max_peak = current_peak;
    }

    // in-set, zero-expressed genes
    for (int j = tid; j < gset_nonzeros; j += GSVA_WALK_THREADS) {
        int key = (nnz - 1) + (int)round(((double) j) * step + 1.0);
        double stat_running = nonzero_total + (double)j * zero_weight;
        gsva_float_t neg_penalty = (gsva_float_t)(key - (nonzero_events + j)) * neg_step;
        double current_valley = stat_running * inv_stat_total - neg_penalty;
        if (current_valley < thread_min_valley) thread_min_valley = current_valley;

        double current_peak =
            (stat_running + zero_weight) * inv_stat_total - neg_penalty;
        if (current_peak > thread_max_peak) thread_max_peak = current_peak;
    }

    PeakValley pv = {(gsva_float_t)thread_max_peak, (gsva_float_t)thread_min_valley};
    PeakValley result = BlockReducePV(tmp.reduce).Reduce(pv, PVReduceOp());

    if (tid == 0) {
        gsva_float_t final_max_peak = result.peak;
        gsva_float_t final_min_valley = result.valley;
        gsva_float_t res = 0.0f;
        switch (score_type) {
            case 0: res = final_max_peak + final_min_valley; break;
            case 1: res = final_max_peak - final_min_valley; break;
            case 2: res = (final_max_peak > GSVA_FABS(final_min_valley)) ? final_max_peak : final_min_valley; break;
            default: res = -1.0f; break;
        }
        out_es[c * S + gset] = res;
    }
}
#endif
