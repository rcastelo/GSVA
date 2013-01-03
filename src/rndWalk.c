#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>

SEXP
rndWalk_R(SEXP gSetIdxR, SEXP geneRankingR, SEXP jR, SEXP RR, SEXP alphaR) {

  int* gSetIdx = INTEGER(gSetIdxR);
  int  k = length(gSetIdxR);
  int* geneRanking = INTEGER(geneRankingR);
  int  p = length(geneRankingR);
  int  j = INTEGER(jR)[0] - 1; /* in C indices are 0-based !!! */
  int* R = INTEGER(RR);
  double alpha = REAL(alphaR)[0];
  SEXP   walkStatR;
  double walkStat;
  double* stepCDFinGeneSet;
  double* stepCDFoutGeneSet;
  double denomStepCDFinGeneSet = 0.0;
  double denomStepCDFoutGeneSet = (double) (p - k);
  int i;

  stepCDFinGeneSet = Calloc(p, double);
  stepCDFoutGeneSet = Calloc(p, double);

  walkStat = 0.0;
  for (i=0; i < k; i++) {
    denomStepCDFinGeneSet += pow(fabs((double) R[j*p+gSetIdx[i]-1]), alpha); /* in C indices are 0-based !!! */
  }

  for (i=0; i < p; i++) {
    int l=0;

    while (l < k && gSetIdx[l] != geneRanking[i])
      l++;

    if (l < k) {
      stepCDFinGeneSet[i] = i == 0 ? pow(fabs((double) R[j*p+geneRanking[i]-1]), alpha) / denomStepCDFinGeneSet :
                                     stepCDFinGeneSet[i-1] + (pow(fabs((double) R[j*p+geneRanking[i]-1]), alpha) / denomStepCDFinGeneSet); /* in C indices are 0-based !!! */
      stepCDFoutGeneSet[i] = stepCDFoutGeneSet[i-1];
    } else {
      stepCDFinGeneSet[i] = stepCDFinGeneSet[i-1];
      stepCDFoutGeneSet[i] = i == 0 ? 1.0 / denomStepCDFoutGeneSet : stepCDFoutGeneSet[i-1] + (1.0 / denomStepCDFoutGeneSet);
    }

    walkStat += stepCDFinGeneSet[i] - stepCDFoutGeneSet[i];
  }

  Free(stepCDFinGeneSet);
  Free(stepCDFoutGeneSet);

  PROTECT(walkStatR = allocVector(REALSXP, 1));
  REAL(walkStatR)[0] = walkStat;
  UNPROTECT(1); /* walkStatR */

  return(walkStatR);
}
