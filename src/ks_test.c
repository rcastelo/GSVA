#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>

/* to add attributes to R objects from C code */
static SEXP
installAttrib(SEXP, SEXP, SEXP);

static SEXP
installAttrib(SEXP vec, SEXP name, SEXP val)
{
  SEXP s, t;

  if (TYPEOF(vec) == CHARSXP)
    error("cannot set attribute on a CHARSXP");
  PROTECT(vec);
  PROTECT(name);
  PROTECT(val);
  for (s = ATTRIB(vec); s != R_NilValue; s = CDR(s)) {
    if (TAG(s) == name) {
      SETCAR(s, val);
      UNPROTECT(3);
      return val;
    }
  }
  s = Rf_allocList(1);
  SETCAR(s, val);
  SET_TAG(s, name);
  if (ATTRIB(vec) == R_NilValue)
    SET_ATTRIB(vec, s);
  else {
    t = nthcdr(ATTRIB(vec), length(ATTRIB(vec)) - 1);
    SETCDR(t, s);
  }
  UNPROTECT(3);
  return val;
}


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
gsva_rnd_walk(int* gsetidx, int k, int* decordstat, double* symrnkstat, int n,
              double tau, double* walkstat, double* walkstatpos,
              double* walkstatneg) {
  int*    gsetrnk;
  double* stepcdfingeneset;
  int*    stepcdfoutgeneset;

  gsetrnk = R_Calloc(k, int);
  for (int i=0; i < k; i++)
    gsetrnk[i] = decordstat[gsetidx[i]-1];

  stepcdfingeneset = R_Calloc(n, double);  /* assuming zeroes are set */
  stepcdfoutgeneset = R_Calloc(n, int);
  for (int i=0; i < n; i++)
    stepcdfoutgeneset[i] = 1;

  for (int i=0; i < k; i++) {
    /* convert 1-based gene indices to 0-based ! */
    if (tau == 1)
      stepcdfingeneset[gsetrnk[i]-1] = symrnkstat[gsetidx[i]-1];
    else
      stepcdfingeneset[gsetrnk[i]-1] = pow(symrnkstat[gsetidx[i]-1], tau);
    stepcdfoutgeneset[gsetrnk[i]-1] = 0;
  }
  R_Free(gsetrnk);

  for (int i=1; i < n; i++) {
    stepcdfingeneset[i] = stepcdfingeneset[i-1] + stepcdfingeneset[i];
    stepcdfoutgeneset[i] = stepcdfoutgeneset[i-1] + stepcdfoutgeneset[i];
  }

  *walkstatpos = *walkstatneg = NA_REAL;
  if (stepcdfingeneset[n-1] > 0 && stepcdfoutgeneset[n-1] > 0) {
    *walkstatpos = *walkstatneg = 0;
    for (int i=0; i < n; i++) {
      double wlkstat = 0;

      if (walkstat != NULL)
        wlkstat = walkstat[i] = ((double) stepcdfingeneset[i]) /
                                ((double) stepcdfingeneset[n-1]) -
                                ((double) stepcdfoutgeneset[i]) /
                                ((double) stepcdfoutgeneset[n-1]);
      else {
        wlkstat = ((double) stepcdfingeneset[i]) / ((double) stepcdfingeneset[n-1]) -
                  ((double) stepcdfoutgeneset[i]) / ((double) stepcdfoutgeneset[n-1]);
      }

      if (wlkstat > *walkstatpos)
        *walkstatpos = wlkstat;
      if (wlkstat < *walkstatneg)
        *walkstatneg = wlkstat;
    }
  }

  R_Free(stepcdfoutgeneset);
  R_Free(stepcdfingeneset);
}

void
gsva_rnd_walk_nas(int* gsetidx, int k, int* decordstat, double* symrnkstat, int n,
                  double tau, int na_use, int minsize, double* walkstat,
                  double* walkstatpos, double* walkstatneg, int* wna) {
  int*    gsetidx_wonas;
  int*    gsetrnk;
  double* stepcdfingeneset;
  int*    stepcdfoutgeneset;
  int     k_notna = 0;

  gsetidx_wonas = R_Calloc(k, int);
  gsetrnk = R_Calloc(k, int);

  for (int i=0; i < k; i++) {
    if (decordstat[gsetidx[i]-1] != NA_INTEGER) { /* na.rm skips NAs */
      gsetidx_wonas[k_notna] = gsetidx[i];
      gsetrnk[k_notna] = decordstat[gsetidx[i]-1];
      k_notna++;
    } else {
      if (na_use < 3) /* everything or all.obs */
        return;
    }
  }

  *walkstatpos = *walkstatneg = NA_REAL;
  if (k_notna >= minsize) { /* na.rm */
    k = k_notna;

    stepcdfingeneset = R_Calloc(n, double);  /* assuming zeroes are set */
    stepcdfoutgeneset = R_Calloc(n, int);
    for (int i=0; i < n; i++)
      stepcdfoutgeneset[i] = 1;

    for (int i=0; i < k; i++) {
      /* convert 1-based gene indices to 0-based ! */
      if (tau == 1)
        stepcdfingeneset[gsetrnk[i]-1] = symrnkstat[gsetidx_wonas[i]-1];
      else
        stepcdfingeneset[gsetrnk[i]-1] = pow(symrnkstat[gsetidx_wonas[i]-1], tau);
      stepcdfoutgeneset[gsetrnk[i]-1] = 0;
    }

    for (int i=1; i < n; i++) {
      stepcdfingeneset[i] = stepcdfingeneset[i-1] + stepcdfingeneset[i];
      stepcdfoutgeneset[i] = stepcdfoutgeneset[i-1] + stepcdfoutgeneset[i];
    }

    if (stepcdfingeneset[n-1] > 0 && stepcdfoutgeneset[n-1] > 0) {
      *walkstatpos = *walkstatneg = 0;
      for (int i=0; i < n; i++) {
        double wlkstat = 0;

        if (walkstat != NULL)
          wlkstat = walkstat[i] = ((double) stepcdfingeneset[i]) /
                                  ((double) stepcdfingeneset[n-1]) -
                                  ((double) stepcdfoutgeneset[i]) /
                                  ((double) stepcdfoutgeneset[n-1]);
        else {
          wlkstat = ((double) stepcdfingeneset[i]) /
                    ((double) stepcdfingeneset[n-1]) -
                    ((double) stepcdfoutgeneset[i]) /
                    ((double) stepcdfoutgeneset[n-1]);
        }

        if (wlkstat > *walkstatpos)
          *walkstatpos = wlkstat;
        if (wlkstat < *walkstatneg)
          *walkstatneg = wlkstat;
      }
    }

    R_Free(stepcdfoutgeneset);
    R_Free(stepcdfingeneset);

  } else
    *wna = 1;

  R_Free(gsetrnk);
  R_Free(gsetidx_wonas);
}

SEXP
gsva_score_genesets_R(SEXP genesetsidxR, SEXP decordstatR, SEXP symrnkstatR,
                      SEXP maxdiffR, SEXP absrnkR, SEXP tauR, SEXP anynaR,
                      SEXP nauseR, SEXP minsizeR) {
  int      m = length(genesetsidxR);
  int      n = length(decordstatR);
  Rboolean maxdiff=asLogical(maxdiffR);
  Rboolean absrnk=asLogical(absrnkR);
  double   tau=REAL(tauR)[0];
  Rboolean anyna=asLogical(anynaR);
  int      nause=INTEGER(nauseR)[0]; /* everything=1, all.obs=2, na.rm=3 */
  int      minsize=INTEGER(minsizeR)[0];
  int*     decordstat;
  double*  symrnkstat;
  SEXP     esR;
  double*  es;
  int      wna=0;
  Rboolean abort=FALSE;

  PROTECT(genesetsidxR);
  PROTECT(decordstatR);
  PROTECT(symrnkstatR);
  PROTECT(esR = allocVector(REALSXP, m));

  decordstat = INTEGER(decordstatR);
  symrnkstat = REAL(symrnkstatR);
  es = REAL(esR);

  for (int i=0; i < m; i++) {
     SEXP    gsetidxR=VECTOR_ELT(genesetsidxR, i);
     int*    gsetidx;
     int     k = length(gsetidxR);
     double  walkstatpos, walkstatneg;

     gsetidx = INTEGER(gsetidxR);
     if (anyna)
       gsva_rnd_walk_nas(gsetidx, k, decordstat, symrnkstat, n, tau,
                         nause, minsize, NULL, &walkstatpos, &walkstatneg,
                         &wna);
     else
       gsva_rnd_walk(gsetidx, k, decordstat, symrnkstat, n, tau,
                     NULL, &walkstatpos, &walkstatneg);

     es[i] = NA_REAL;
     if (!anyna || (!ISNA(walkstatpos) && !ISNA(walkstatneg))) {
	     if (maxdiff) {
		     es[i] = walkstatpos + walkstatneg;
         if (absrnk)
           es[i] = walkstatpos - walkstatneg;
	     } else {
		       es[i] = (walkstatpos > fabs(walkstatneg)) ? walkstatpos : walkstatneg;
	     }
     } else {
       if (anyna && (ISNA(walkstatpos) || ISNA(walkstatneg)) && nause == 2) { /* all.obs */
         abort=TRUE;
         break;
       }
     }
  }

  if (anyna) {
    SEXP class;

    if (nause == 2 && abort) {
      PROTECT(class = allocVector(STRSXP, 1));
      SET_STRING_ELT(class, 0, mkChar("abort"));
      installAttrib(esR, R_ClassSymbol, class);
      UNPROTECT(1); /* class */
    } else if (nause == 3 && wna == 1) {
      PROTECT(class = allocVector(STRSXP, 1));
      SET_STRING_ELT(class, 0, mkChar("wna"));
      installAttrib(esR, R_ClassSymbol, class);
      UNPROTECT(1); /* class */
    }
  }

  UNPROTECT(4); /* genesetsidxR decordstatR symrnkstatR esR */

  return(esR);
}
