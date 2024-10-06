#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>

SEXP
ks_matrix_R(SEXP XR, SEXP sidxsR, SEXP n_genesR, SEXP geneset_idxsR,
            SEXP n_genesetR, SEXP tauR, SEXP n_samplesR, SEXP mx_diffR, SEXP abs_rnkR);

double
ks_sample(double* x, int* x_sort_indxs, int n_genes, int* geneset_mask,
          int* geneset_idxs, int n_geneset, double tau, int mx_diff, int abs_rnk){


	double dec = 1.0 / (n_genes - n_geneset);
	double sum_gset = 0.0;
	for(int i = 0; i < n_geneset; ++i){
		sum_gset += pow(x[geneset_idxs[i]-1], tau);
	}
	//Rprintf("%d\n", mx_diff);

	//double mx_value = 0.0;
	double mx_value_sign = 0.0;
	double cum_sum = 0.0;

	double mx_pos = 0.0;
	double mx_neg = 0.0;

	int idx;
	for(int i = 0; i < n_genes; ++i){
		idx = x_sort_indxs[i]-1;

		if(geneset_mask[idx] == 1){
			cum_sum += pow(x[idx], tau) / sum_gset;
		}else{
			cum_sum -= dec;
		}

		if(cum_sum > mx_pos){ mx_pos = cum_sum; }
		if(cum_sum < mx_neg){ mx_neg = cum_sum; }
	}

	if (mx_diff != 0) {
		mx_value_sign = mx_pos + mx_neg;
    if (abs_rnk != 0)
      mx_value_sign = mx_pos - mx_neg;
	} else {
		mx_value_sign = (mx_pos > fabs(mx_neg)) ? mx_pos : mx_neg;
	}
	return mx_value_sign;
}


/**
 * X <- gene density scores
 * R <- result
 * sidxs <- sorted gene densities idxs
 */
void ks_matrix(double* X, double* R, int* sidxs, int n_genes, int* geneset_idxs,
               int n_geneset, double tau, int n_samples, int mx_diff, int abs_rnk){
	int geneset_mask[n_genes];
	for(int i = 0; i < n_genes; ++i){
		geneset_mask[i] = 0;
	}

	for(int i = 0; i < n_geneset; ++i){
		geneset_mask[geneset_idxs[i]-1] = 1;
	}

	for(int j = 0; j < n_samples; ++j){
		int offset = j * n_genes;
    R[j] = ks_sample(&X[offset], &sidxs[offset], n_genes, &geneset_mask[0],
                     geneset_idxs, n_geneset, tau, mx_diff, abs_rnk);
	}
}

SEXP
ks_matrix_R(SEXP XR, SEXP sidxsR, SEXP n_genesR, SEXP geneset_idxsR,
            SEXP n_genesetR, SEXP tauR, SEXP n_samplesR, SEXP mx_diffR,
            SEXP abs_rnkR) {
  double* X=REAL(XR);
  int*    sidxs=INTEGER(sidxsR);
  int     n_genes=INTEGER(n_genesR)[0];
  int*    geneset_idxs=INTEGER(geneset_idxsR);
  int     n_geneset=INTEGER(n_genesetR)[0];
  double  tau=REAL(tauR)[0];
  int     n_samples=INTEGER(n_samplesR)[0];
  int     mx_diff=INTEGER(mx_diffR)[0];
  int     abs_rnk=INTEGER(abs_rnkR)[0];  
  SEXP    resR;
  double* res;

  PROTECT(resR = allocVector(REALSXP, n_samples));
  res = REAL(resR);

  ks_matrix(X, res, sidxs, n_genes, geneset_idxs, n_geneset, tau, n_samples,
            mx_diff, abs_rnk);

  UNPROTECT(1); /* resR */

  return(resR);
}

void
gsva_rnd_walk(int* gsetidx, int k, int* generanking, int* rankstat, int n,
              double* walkstat, double* walkstatpos, double* walkstatneg) {
  int*    stepcdfingeneset;
  int*    stepcdfoutgeneset;

  stepcdfingeneset = R_Calloc(n, int);  /* assuming zeroes are set */
  stepcdfoutgeneset = R_Calloc(n, int);
  for (int i=0; i < n; i++)
    stepcdfoutgeneset[i] = 1;

  for (int i=0; i < k; i++) {
    /* convert 1-based gene indices to 0-based ! */
    stepcdfingeneset[gsetidx[i]-1] = rankstat[generanking[gsetidx[i]-1]-1];
    stepcdfoutgeneset[gsetidx[i]-1] = 0;
  }

  for (int i=1; i < n; i++) {
    stepcdfingeneset[i] = stepcdfingeneset[i-1] + stepcdfingeneset[i];
    stepcdfoutgeneset[i] = stepcdfoutgeneset[i-1] + stepcdfoutgeneset[i];
  }

  *walkstatpos = *walkstatneg = 0;
  for (int i=0; i < n; i++) {
    double wlkstat = 0;

    if (walkstat != NULL)
      wlkstat = walkstat[i] = ((double) stepcdfingeneset[i]) / ((double) stepcdfingeneset[n-1]) -
                              ((double) stepcdfoutgeneset[i]) / ((double) stepcdfoutgeneset[n-1]);
    else {
      wlkstat = ((double) stepcdfingeneset[i]) / ((double) stepcdfingeneset[n-1]) -
                ((double) stepcdfoutgeneset[i]) / ((double) stepcdfoutgeneset[n-1]);
    }

    if (wlkstat > *walkstatpos)
      *walkstatpos = wlkstat;
    if (wlkstat < *walkstatneg)
      *walkstatneg = wlkstat;
  }

  R_Free(stepcdfoutgeneset);
  R_Free(stepcdfingeneset);
}

SEXP
gsva_rnd_walk_R(SEXP gsetidxR, SEXP generankingR, SEXP rankstatR) {
  int     n = length(generankingR);
  int     k = length(gsetidxR);
  int*    gsetidx;
  int*    generanking;
  int*    rankstat;
  SEXP    walkstatR;
  double* walkstat;
  double  walkstatpos, walkstatneg;

  PROTECT(gsetidxR);
  PROTECT(generankingR);
  PROTECT(rankstatR);
  PROTECT(walkstatR = allocVector(REALSXP, n));

  gsetidx = INTEGER(gsetidxR);
  generanking = INTEGER(generankingR);
  rankstat = INTEGER(rankstatR);
  walkstat = REAL(walkstatR);

  gsva_rnd_walk(gsetidx, k, generanking, rankstat, n,
                walkstat, &walkstatpos, &walkstatneg);

  UNPROTECT(4); /* gsetidxR generankingR rankstatR walkstatR */

  return(walkstatR);
}

void
gsva_rnd_walk_nonunittau(int* gsetidx, int k, int* generanking, int* rankstat,
                         int n, double tau,
                         double* walkstat, double* walkstatpos, double* walkstatneg) {
  double* stepcdfingeneset;
  double* stepcdfoutgeneset;

  stepcdfingeneset = R_Calloc(n, double);  /* assuming zeroes are set */
  stepcdfoutgeneset = R_Calloc(n, double);
  for (int i=0; i < n; i++)
    stepcdfoutgeneset[i] = 1.0;

  for (int i=0; i < k; i++) {
    /* convert 1-based gene indices to 0-based ! */
    stepcdfingeneset[gsetidx[i]-1] = pow(rankstat[generanking[gsetidx[i]-1]-1], tau);
    stepcdfoutgeneset[gsetidx[i]-1] = 0;
  }

  for (int i=1; i < n; i++) {
    stepcdfingeneset[i] = stepcdfingeneset[i-1] + stepcdfingeneset[i];
    stepcdfoutgeneset[i] = stepcdfoutgeneset[i-1] + stepcdfoutgeneset[i];
  }

  *walkstatpos = *walkstatneg = 0;
  for (int i=0; i < n; i++) {
    double wlkstat = 0;

    if (walkstat != NULL)
      wlkstat = walkstat[i] = ((double) stepcdfingeneset[i]) / ((double) stepcdfingeneset[n-1]) -
                              ((double) stepcdfoutgeneset[i]) / ((double) stepcdfoutgeneset[n-1]);
    else {
      wlkstat = ((double) stepcdfingeneset[i]) / ((double) stepcdfingeneset[n-1]) -
                ((double) stepcdfoutgeneset[i]) / ((double) stepcdfoutgeneset[n-1]);
    }

    if (wlkstat > *walkstatpos)
      *walkstatpos = wlkstat;
    if (wlkstat < *walkstatneg)
      *walkstatneg = wlkstat;
  }

  R_Free(stepcdfoutgeneset);
  R_Free(stepcdfingeneset);
}

SEXP
gsva_score_genesets_R(SEXP genesetsrankidxR, SEXP generankingR, SEXP rankstatR,
                      SEXP maxdiffR, SEXP absrnkR, SEXP tauR) {
  int      m = length(genesetsrankidxR);
  int      n = length(generankingR);
  Rboolean maxdiff=asLogical(maxdiffR);
  Rboolean absrnk=asLogical(absrnkR);
  double   tau=REAL(tauR)[0];
  int*     generanking;
  int*     rankstat;
  SEXP     esR;
  double*  es;

  PROTECT(genesetsrankidxR);
  PROTECT(generankingR);
  PROTECT(rankstatR);
  PROTECT(esR = allocVector(REALSXP, m));

  generanking = INTEGER(generankingR);
  rankstat = INTEGER(rankstatR);
  es = REAL(esR);

  for (int i=0; i < m; i++) {
     SEXP    gsetidxR=VECTOR_ELT(genesetsrankidxR, i);
     int*    gsetidx;
     int     k = length(gsetidxR);
     double  walkstatpos, walkstatneg;

     gsetidx = INTEGER(gsetidxR);
     if (tau == 1) /* treated separately to reduce memory consumption */
       gsva_rnd_walk(gsetidx, k, generanking, rankstat, n,
                     NULL, &walkstatpos, &walkstatneg);
     else
       gsva_rnd_walk_nonunittau(gsetidx, k, generanking, rankstat, n, tau,
                     NULL, &walkstatpos, &walkstatneg);

	   if (maxdiff) {
		   es[i] = walkstatpos + walkstatneg;
       if (absrnk)
         es[i] = walkstatpos - walkstatneg;
	   } else {
		     es[i] = (walkstatpos > fabs(walkstatneg)) ? walkstatpos : walkstatneg;
	   }
  }

  UNPROTECT(4); /* genesetsrankidxR generankingR rankstatR esR */

  return(esR);
}
