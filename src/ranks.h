#include <R.h>
#include <Rdefines.h>

/* type codes */
#define RANKSTYPE_MATRIX_INT   1
#define RANKSTYPE_MATRIX_DBL   2
#define RANKSTYPE_DGC          3
#define RANKSTYPE_SVT_INT      4
#define RANKSTYPE_SVT_DBL      5

typedef struct ranks_ctx_s {
  int type;            /* RANKSTYPE_*  */
  union {
    const void* data;  /* matrix: INTEGER or REAL base ptr     */
    struct { const int* i; const int* p; const double* x; } dgc;
    SEXP svt;
  } u;
  int p;
  int n;
  Rboolean sparse;
  int (*fetch_col)(struct ranks_ctx_s* ctx, int j, int* col);
} ranks_ctx_t;

ranks_ctx_t* ranks_ctx_create(SEXP XR, Rboolean intrnks, Rboolean sparse);

int
ranks2stats(ranks_ctx_t* ctx, int j, int* decordstat, double* symrnkstat,
            int* nzeros_col, double* zerosymrnkstat);

int
ranks2stats_nas(ranks_ctx_t* ctx, int j, int* decordstat, double* symrnkstat,
                int* nzeros_col, double* zerosymrnkstat, int* nnas);
