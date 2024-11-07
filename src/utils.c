#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <cli/progress.h>

/* global variables */
extern SEXP Matrix_DimNamesSym,
            Matrix_DimSym,
            Matrix_xSym,
            Matrix_iSym,
            Matrix_jSym,
            Matrix_pSym;

/* indirect numerical (double) comparison for qsort() for obtaining
 * permutations leading to a decreasing order
 */
double* global_dbl_p;

double
sd(double* x, int n);

double
sd_naprop(double* x, int n);

double
sd_narm(double* x, int n);


/* calculates standard deviation, largely borrowed from C code in R's src/main/cov.c */
double
sd(double* x, int n) {
  int         i, n1;
  double      mean, sd;
  long double sum = 0.0;
  long double tmp;

  for (i=0; i < n; i++)
    sum += x[i];
  tmp = sum / n;
  if (R_FINITE((double) tmp)) {
    sum = 0.0;
    for (i=0; i < n; i++)
      sum += x[i] - tmp;
    tmp = tmp + sum / n;
  }
  mean = tmp;
  n1 = n - 1;

  sum = 0.0;
  for (i=0; i < n; i++)
    sum += (x[i] - mean) * (x[i] - mean);
  sd = sqrt((double) (sum / ((long double) n1)));

  return(sd);
}

/* calculates standard deviation, largely borrowed from C code in R's
 * src/main/cov.c propagating NA values */
double
sd_naprop(double* x, int n) {
  int         i, n1;
  double      mean, sd;
  long double sum = 0.0;
  long double tmp;

  for (i=0; i < n; i++) {
    if (!ISNA(x[i]))
      sum += x[i];
    else
      return(NA_REAL);
  }
  tmp = sum / n;
  if (R_FINITE((double) tmp)) {
    sum = 0.0;
    for (i=0; i < n; i++)
      sum += x[i] - tmp;
    tmp = tmp + sum / n;
  }
  mean = tmp;
  n1 = n - 1;

  sum = 0.0;
  for (i=0; i < n; i++)
    sum += (x[i] - mean) * (x[i] - mean);
  sd = sqrt((double) (sum / ((long double) n1)));

  return(sd);
}

/* calculates standard deviation, largely borrowed from C code in R's
 * src/main/cov.c removing NA values */
double
sd_narm(double* x, int n) {
  int         i, n1;
  double      mean, sd;
  long double sum = 0.0;
  long double tmp;
  int         n_nas = 0;

  for (i=0; i < n; i++) {
    if (!ISNA(x[i]))
      sum += x[i];
    else
      n_nas++;
  }
  if (n_nas >= n - 1)
    return(NA_REAL);

  tmp = sum / (n - n_nas);
  if (R_FINITE((double) tmp)) {
    sum = 0.0;
    for (i=0; i < n; i++)
      if (!ISNA(x[i]))
        sum += x[i] - tmp;
    tmp = tmp + sum / (n - n_nas);
  }
  mean = tmp;
  n1 = n - n_nas - 1;

  sum = 0.0;
  for (i=0; i < n; i++)
    if (!ISNA(x[i]))
      sum += (x[i] - mean) * (x[i] - mean);
  sd = sqrt((double) (sum / ((long double) n1)));

  return(sd);
}

int
indirect_dbl_cmp_dec(const void* a, const void* b) {
  const int *ai = (int *) a;
  const int *bi = (int *) b;

  /* assuming index are 1-based! */
  if (global_dbl_p[*ai-1] < global_dbl_p[*bi-1])
    return 1;
  if (global_dbl_p[*ai-1] > global_dbl_p[*bi-1])
    return -1;

  return 0;
}

/* src/main/unique.c */
/*
SEXP match5(SEXP itable, SEXP ix, int nmatch, SEXP incomp, SEXP env);
*/

SEXP
match_int(SEXP x, SEXP table);

SEXP
match_int(SEXP x, SEXP table) {
  SEXP s, t, ansR;

  PROTECT(t = s = LCONS(R_NilValue, Rf_allocList(3)));
  SETCAR(t, install("match")); t=CDR(t);
  SETCAR(t, x);
  SET_TAG(t, install("x")); t=CDR(t);
  SETCAR(t, table);
  SET_TAG(t, install("table"));
  ansR = eval(s, R_GlobalEnv);

  /*
  SEXP ansR = match5(table, x, NA_INTEGER, NULL, R_GlobalEnv);
  */
  UNPROTECT(1); /* t s */

  return(ansR);
}

SEXP
order_rankstat_R(SEXP xR);

void
order_rankstat(double* x, int n, int* ord, double* rst) {
  for (int i=0; i < n; i++)
    ord[i] = i + 1; /* 1-based to comply with R upstream */

  global_dbl_p = x;
  qsort(ord, n, sizeof(int), indirect_dbl_cmp_dec);

  for (int i=0; i < n; i++)
    rst[ord[i]-1] = fabs(((double) n) - ((double) i) - (((double ) n) / 2.0));
}

SEXP
order_rankstat_R(SEXP xR) {
  int  n = length(xR);
  double* x;
  SEXP ordR, rstR, ansR;
  int* ord;
  double* rst;

  PROTECT(xR);
  x = REAL(xR);

  PROTECT(ordR = allocVector(INTSXP, n));
  PROTECT(rstR = allocVector(REALSXP, n));

  ord = INTEGER(ordR);
  rst = REAL(rstR);

  order_rankstat(x, n, ord, rst);

  PROTECT(ansR = allocVector(VECSXP, 2));

  SET_VECTOR_ELT(ansR, 0, ordR);
  SET_VECTOR_ELT(ansR, 1, rstR);

  UNPROTECT(4); /* xR ordR rstR ansR */

  return(ansR);
}

SEXP
order_rankstat_sparse_to_dense_R(SEXP XCspR, SEXP jR) { /* column in jR is 1-based! */
  int  j = INTEGER(jR)[0];
  int* XCsp_dim;
  int* XCsp_i;
  int* XCsp_p;
  double* XCsp_x;
  double* x;
  int     nr;
  SEXP    ordR, rstR, ansR;
  int*    ord;
  double* rst;

  PROTECT(XCspR);
  XCsp_dim = INTEGER(GET_SLOT(XCspR, Matrix_DimSym));
  nr = XCsp_dim[0]; /* number of rows */
  XCsp_i = INTEGER(GET_SLOT(XCspR, Matrix_iSym));
  XCsp_p = INTEGER(GET_SLOT(XCspR, Matrix_pSym));
  XCsp_x = REAL(GET_SLOT(XCspR, Matrix_xSym));

  /* put the sparse column into a dense vector */
  x = R_Calloc(nr, double);
  for (int i=XCsp_p[j-1]; i < XCsp_p[j]; i++)
    x[XCsp_i[i]] = XCsp_x[i];

  PROTECT(ordR = allocVector(INTSXP, nr));
  PROTECT(rstR = allocVector(REALSXP, nr));

  ord = INTEGER(ordR);
  rst = REAL(rstR);

  order_rankstat(x, nr, ord, rst);

  R_Free(x);

  PROTECT(ansR = allocVector(VECSXP, 2));

  SET_VECTOR_ELT(ansR, 0, ordR);
  SET_VECTOR_ELT(ansR, 1, rstR);

  UNPROTECT(4); /* XCspR ordR rstR ansR */

  return(ansR);
}

/* here the 'sparse_to_sparse' suffix does not refer to returning
 * a sparse data structure, but to apply a sparse strategy, the
 * returned value is anyway dense */
SEXP
order_rankstat_sparse_to_sparse_R(SEXP XCspR, SEXP jR) { /* column in jR is 1-based! */
  int  j = INTEGER(jR)[0];
  int* XCsp_dim;
  int* XCsp_i;
  int* XCsp_p;
  double* XCsp_x;
  int     nnz_j;
  double* x;
  int     nr;
  SEXP    ordR, rstR, ansR;
  int*    ord;
  double* rst;
  int*    ord2;
  int*    allord;
  int     allord_i;

  PROTECT(XCspR);
  XCsp_dim = INTEGER(GET_SLOT(XCspR, Matrix_DimSym));
  nr = XCsp_dim[0]; /* number of rows */
  XCsp_i = INTEGER(GET_SLOT(XCspR, Matrix_iSym));
  XCsp_p = INTEGER(GET_SLOT(XCspR, Matrix_pSym));
  XCsp_x = REAL(GET_SLOT(XCspR, Matrix_xSym));

  /* put the nonzero values of the sparse column into a smaller dense vector */
  nnz_j = XCsp_p[j] - XCsp_p[j-1];
  x = R_Calloc(nnz_j, double);
  for (int i=XCsp_p[j-1]; i < XCsp_p[j]; i++) {
    int k = i - XCsp_p[j-1];
    x[k] = XCsp_x[i];
  }

  allord = R_Calloc(nr, int);
  for (int i=0; i < nr; i++)
    allord[i] = i + 1;

  PROTECT(ordR = allocVector(INTSXP, nr));
  PROTECT(rstR = allocVector(REALSXP, nr));

  ord2 = R_Calloc(nnz_j, int);
  ord = INTEGER(ordR);
  rst = REAL(rstR);

  for (int i=0; i < nnz_j; i++)
    ord2[i] = i+1; /* indirect_dbl_cmp_dec() assumes 1-based! */

  global_dbl_p = x;
  qsort(ord2, nnz_j, sizeof(int), indirect_dbl_cmp_dec);

  for (int i=0; i < nnz_j; i++) {
    ord[i] = XCsp_i[ord2[i] - 1 + XCsp_p[j-1]] + 1; /* should return 1-based ! */
    allord[ord[i]-1] = -1;
  }

  allord_i = nnz_j;
  for (int i=0; i < nr; i++) /* place zeroes at the end of the ranking */
    if (allord[i] > 0) {
      ord[allord_i] = allord[i];
      allord_i++;
    }

  for (int i=0; i < nr; i++)
    rst[i] = (nnz_j / 2.0) + 1.0; /* zero entries get the same symmetric rank */

  for (int i=0; i < nnz_j; i++)
    rst[ord[i]-1] = fabs(((double) nnz_j) - ((double) i) - (((double) nnz_j) / 2.0));

  R_Free(ord2);
  R_Free(allord);
  R_Free(x);

  PROTECT(ansR = allocVector(VECSXP, 2));

  SET_VECTOR_ELT(ansR, 0, ordR);
  SET_VECTOR_ELT(ansR, 1, rstR);

  UNPROTECT(4); /* XCspR ordR rstR ansR */

  return(ansR);
}


/* from https://github.com/cran/curl/blob/master/src/interrupt.c */
/* Check for interrupt without long jumping */
void
check_interrupt_fn(void *dummy) {
  R_CheckUserInterrupt();
}

int
pending_interrupt(void) {
  return !(R_ToplevelExec(check_interrupt_fn, NULL));
}
