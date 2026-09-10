#include <R.h>
#include <Rdefines.h>

typedef int (*FetchColFunDef)(SEXP, int, int, int*);

FetchColFunDef find_dim_and_fetchcolfun(SEXP XR, Rboolean intrnks, int** dim);

int
ranks2stats(SEXP ranksR, int p, int n, int j, Rboolean sparse,
            FetchColFunDef fetch_col,
            int* decordstat_col, double* symrnkstat_col,
            int* nzeros_col, double* zerosymrnkstat);

int
ranks2stats_nas(SEXP ranksR, int p, int n, int j, Rboolean sparse,
                FetchColFunDef fetch_col,
                int* decordstat_col, double* symrnkstat_col,
                int* nzeros_col, double* zerosymrnkstat, int* nnas);
