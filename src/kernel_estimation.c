/*
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
*/
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <cli/progress.h>

SEXP
matrix_density_R(SEXP density_dataR, SEXP test_dataR, SEXP n_density_samplesR,
                 SEXP n_test_samplesR, SEXP n_genesR, SEXP GausskR,
                 SEXP any_naR, SEXP na_useR, SEXP verboseR);

void initCdfs(void);
double precomputedCdf(double x, double sigma);

#define SIGMA_FACTOR 4.0
#define PRECOMPUTE_RESOLUTION 10000
#define MAX_PRECOMPUTE 10.0

double precomputed_cdf[PRECOMPUTE_RESOLUTION+1];
int is_precomputed = 0;

double sd(double* x, int n);

double sd_naprop(double* x, int n);

double sd_narm(double* x, int n);

/* including expression log-odds */
void
row_d(double* x, double* y, double* r, int size_density_n,
      int size_test_n, int Gaussk) {

  double bw = Gaussk ? (sd(x, size_density_n) / SIGMA_FACTOR) : 0.5;

  if (ISNA(bw) || bw == 0) /* just in case */
    bw = 0.001;

  if (Gaussk && is_precomputed == 0) {
    initCdfs();
    is_precomputed = 1;
  }

	for(int j = 0; j < size_test_n; ++j){
		double left_tail = 0.0;

		for(int i = 0; i < size_density_n; ++i){
			left_tail += Gaussk ? precomputedCdf(y[j]-x[i], bw) : ppois(y[j], x[i]+bw, TRUE, FALSE);
		}
		left_tail = left_tail / size_density_n;
		r[j] = -1.0 * log((1.0-left_tail)/left_tail);
	}
}

/* including expression log-odds, propagating NAs */
void
row_d_naprop(double* x, double* y, double* r, int size_density_n,
             int size_test_n, int Gaussk) {

  double bw = Gaussk ? (sd_naprop(x, size_density_n) / SIGMA_FACTOR) : 0.5;

  if (!ISNA(bw) && bw == 0) /* just in case */
    bw = 0.001;

  if (Gaussk && is_precomputed == 0) {
    initCdfs();
    is_precomputed = 1;
  }

	for (int j = 0; j < size_test_n; ++j) {
		double left_tail = 0.0;

    if (!ISNA(bw) && !ISNA(y[j])) {
      int i = 0;
      int noISNA = 1;
      while (noISNA && i < size_density_n) {
          noISNA = !ISNA(x[i]);
          if (noISNA) {
			      left_tail += Gaussk ? precomputedCdf(y[j]-x[i], bw) :
                                  ppois(y[j], x[i]+bw, TRUE, FALSE);
            i++;
          }
		  }
      r[j] = NA_REAL;
      if (noISNA) {
		    left_tail = left_tail / size_density_n;
		    r[j] = -1.0 * log((1.0-left_tail)/left_tail);
      }
    } else
      r[j] = NA_REAL;
	}
}

void
row_d_narm(double* x, double* y, double* r, int size_density_n,
           int size_test_n, int Gaussk) {

  double bw = Gaussk ? (sd_narm(x, size_density_n) / SIGMA_FACTOR) : 0.5;

  if (!ISNA(bw) && bw == 0) /* just in case */
    bw = 0.001;

  if (Gaussk && is_precomputed == 0) {
    initCdfs();
    is_precomputed = 1;
  }

	for(int j = 0; j < size_test_n; ++j){
		double left_tail = 0.0;
    int    n_nas = 0;

    if (!ISNA(bw) && !ISNA(y[j])) {
		  for(int i = 0; i < size_density_n; ++i){
        if (!ISNA(x[i]))
			    left_tail += Gaussk ? precomputedCdf(y[j]-x[i], bw) : ppois(y[j], x[i]+bw, TRUE, FALSE);
        else
          n_nas++;
		  }
      if (n_nas < size_density_n) {
		    left_tail = left_tail / (size_density_n - n_nas);
		    r[j] = -1.0 * log((1.0-left_tail)/left_tail);
      } else
        r[j] = NA_REAL;
    } else
      r[j] = NA_REAL;
	}
}

/* without expression log-odds, called from sparse methods */
void
row_d_nologodds(double* x, double* y, double* r, int size_density_n,
                int size_test_n, int Gaussk);
void
row_d_nologodds(double* x, double* y, double* r, int size_density_n,
                int size_test_n, int Gaussk) {

  double bw = Gaussk ? (sd(x, size_density_n) / SIGMA_FACTOR) : 0.5;

  if (ISNA(bw) || bw == 0) /* just in case */
    bw = 0.001;

  if (Gaussk && is_precomputed == 0) {
    initCdfs();
    is_precomputed = 1;
  }

	for (int j = 0; j < size_test_n; ++j) {
		double left_tail = 0.0;

		for (int i = 0; i < size_density_n; ++i) {
      if (ISNA(x[i]) || ISNA(y[j]) || ISNAN(x[i]) || ISNAN(y[j]))
          error("Sparse nonzero expression values cannot be NAs or NaNs.");
			left_tail += Gaussk ? precomputedCdf(y[j]-x[i], bw) : ppois(y[j], x[i]+bw, TRUE, FALSE);
		}
		r[j] = left_tail / size_density_n;
	}
}


void
matrix_d(double* X, double* Y, double* R, int n_density_samples,
         int n_test_samples, int n_genes, int Gaussk, Rboolean any_na,
         int na_use, Rboolean verbose) {
  SEXP pb = R_NilValue;

  if (verbose) {
    pb = PROTECT(cli_progress_bar(n_genes, NULL));
    cli_progress_set_name(pb, "Estimating ECDFs");
  }
    
	for(int j = 0; j < n_genes; ++j){
#ifdef LONG_VECTOR_SUPPORT
    R_xlen_t offset_density = (R_xlen_t) n_density_samples * j;
    R_xlen_t offset_test = (R_xlen_t) n_test_samples * j;
#else
		int      offset_density = (size_t) n_density_samples * j;
		int      offset_test = (size_t) n_test_samples * j;
#endif

    if (!any_na)
		  row_d(&X[offset_density], &Y[offset_test], &R[offset_test],
            n_density_samples, n_test_samples, Gaussk);
    else {
      if (na_use == 1L) /* propagate NAs */
		    row_d_naprop(&X[offset_density], &Y[offset_test], &R[offset_test],
                     n_density_samples, n_test_samples, Gaussk);
      else              /* remove NAs (assuming 3 b/c 2 should have prompt
                           the error earlier */
		    row_d_narm(&X[offset_density], &Y[offset_test], &R[offset_test],
                   n_density_samples, n_test_samples, Gaussk);
    }

    if (verbose) { /* show progress */
      if (j % 100 == 0 && CLI_SHOULD_TICK)
        cli_progress_set(pb, j);
    }
	}

  if (verbose) {
    cli_progress_done(pb);
    UNPROTECT(1); /* pb */
  }
}

SEXP
matrix_density_R(SEXP density_dataR, SEXP test_dataR, SEXP n_density_samplesR,
                 SEXP n_test_samplesR, SEXP n_genesR, SEXP GausskR,
                 SEXP any_naR, SEXP na_useR, SEXP verboseR) {
  double* density_data=REAL(density_dataR);
  double* test_data=REAL(test_dataR);
  int     n_density_samples=INTEGER(n_density_samplesR)[0];
  int     n_test_samples=INTEGER(n_test_samplesR)[0];
  int     n_genes=INTEGER(n_genesR)[0];
  int     Gaussk=INTEGER(GausskR)[0];
  Rboolean any_na=asLogical(any_naR);
  int      na_use=INTEGER(na_useR)[0]; /* everything=1 all.obs=2 na.rm=3 */
  Rboolean verbose=asLogical(verboseR);
  SEXP    resR;
  double* res;

  PROTECT(resR = allocVector(REALSXP, n_test_samples * n_genes));
  res = REAL(resR);

  matrix_d(density_data, test_data, res, n_density_samples, n_test_samples,
           n_genes, Gaussk, any_na, na_use, verbose);

  UNPROTECT(1); /* resR */

  return(resR);
}

inline double precomputedCdf(double x, double sigma){
	double v = x / sigma;
	if(v < (-1 * MAX_PRECOMPUTE)){
		return 0;
	}else if(v > MAX_PRECOMPUTE){
		return 1;
	}else{
		double cdf = precomputed_cdf[(int)(fabs(v) / MAX_PRECOMPUTE * PRECOMPUTE_RESOLUTION)];
		if(v < 0){
			return 1.0 - cdf;
		}else{
			return cdf;
		}
	}
}

void initCdfs(void){
	double divisor = PRECOMPUTE_RESOLUTION * 1.0;
	for(int i = 0; i <= PRECOMPUTE_RESOLUTION; ++i)
    precomputed_cdf[i] = pnorm5(MAX_PRECOMPUTE * ((double) i) / divisor, 0.0, 1.0, TRUE, FALSE);
                         /* standard normal distribution function, lower.tail=TRUE, log.p=FALSE */
}
