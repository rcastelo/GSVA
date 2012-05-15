/*
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
*/
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>

/**
 * X <- possibly bootstrapped samples for density distribution
 * Y <- samples to test
 */
void
matrix_density_R(double* X, double* Y, double* R, int* n_density_samples, int* n_test_samples, int* n_genes, int* rnaseq);

void initCdfs(void);
double precomputedCdf(double x, double sigma);

#define SIGMA_FACTOR 4.0
#define PRECOMPUTE_RESOLUTION 10000
#define MAX_PRECOMPUTE 10.0

double precomputed_cdf[PRECOMPUTE_RESOLUTION+1];
int is_precomputed = 0;

/* calculates standard deviation, largely borrowed from C code in R's src/main/cov.c */
double sd(double* x, int n) {
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

/**
 * for resampling, x are the resampled points and y are the
 */
void row_d(double* x, double* y, double* r, int size_density_n, int size_test_n, int rnaseq){

  double bw = rnaseq ? 0.5 : (sd(x, size_density_n) / SIGMA_FACTOR);

  if (!rnaseq && is_precomputed == 0) {
    initCdfs();
    is_precomputed = 1;
  }

	for(int j = 0; j < size_test_n; ++j){
		double left_tail = 0.0;

		for(int i = 0; i < size_density_n; ++i){
			left_tail += rnaseq ? ppois(y[j], x[i]+bw, TRUE, FALSE) : precomputedCdf(y[j]-x[i], bw);
		}
		left_tail = left_tail / size_density_n;
		/* r[j] = -1.0 * log((1.0-left_tail)/left_tail); */ /* log-odds not necessary */
    r[j] = left_tail;
	}
}

void matrix_d(double* X, double* Y, double* R, int n_density_samples, int n_test_samples, int n_genes, int rnaseq){
	for(int j = 0; j < n_genes; ++j){
		int offset_density = j * n_density_samples;
		int offset_test = j * n_test_samples;
		row_d(&X[offset_density], &Y[offset_test], &R[offset_test], n_density_samples, n_test_samples, rnaseq);
	}
}

void matrix_density_R(double* density_data, double* test_data, double* R, int* n_density_samples,
                      int* n_test_samples, int* n_genes, int* rnaseq){
	matrix_d(density_data, test_data, R,*n_density_samples,*n_test_samples, *n_genes, *rnaseq);
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


