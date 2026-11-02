#ifndef RANKS_H
#define RANKS_H

#include <R.h>
#include <Rdefines.h>

/* type codes of the rank matrix */
#define RANKSTYPE_MATRIX_INT   1
#define RANKSTYPE_MATRIX_DBL   2
#define RANKSTYPE_DGC          3
#define RANKSTYPE_SVT_INT      4
#define RANKSTYPE_SVT_DBL      5

/* context for fetching the columns of a rank matrix, which is created once
 * per call to gsva_score_genesets_R() and holds a pointer to the data of the
 * rank matrix, according to its class, and a buffer of 'p' integers that is
 * reused to fetch every column */
typedef struct ranks_ctx_s {
  int      type;       /* RANKSTYPE_* */
  union {
    const void* data;  /* dense matrix: INTEGER() or REAL() pointer */
    struct { const int* i; const int* p; const double* x; } dgc;
    SEXP svt;          /* SVT slot of a SVT_SparseMatrix */
  } u;
  int      p;          /* number of rows/genes/features */
  int      n;          /* number of columns/samples/cells/spots */
  Rboolean sparse;
  int    (*fetch_col)(struct ranks_ctx_s* ctx, int j);
  int*     r;          /* buffer with the last fetched column */
} ranks_ctx_t;

const char*
get_class_name(SEXP obj);

ranks_ctx_t*
ranks_ctx_create(SEXP XR, Rboolean intrnks, Rboolean sparse);

int
ranks2stats(ranks_ctx_t* ctx, int j, int* decordstat, double* symrnkstat,
            int* nzeros_col, double* zerosymrnkstat);

int
ranks2stats_nas(ranks_ctx_t* ctx, int j, int* decordstat, double* symrnkstat,
                int* nzeros_col, double* zerosymrnkstat, int* nnas);

#endif
