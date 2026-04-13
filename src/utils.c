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

  if (n < 2)
    return(NA_REAL);

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

  if (n < 2)
    return(NA_REAL);

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
