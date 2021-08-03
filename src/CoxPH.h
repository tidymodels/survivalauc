#ifndef _COXPH_H_
#define _COXPH_H_

void DataM(int n, double *Y, double *D, int *L, double *uniqY, int *locY);

void coxph(int n, double *Y, double *Delta, int p, double **X, double *beta, double *lambda,
           double *Lambda, int *L, double *uniqY, int *locY, double *betase, double **cov);

#endif
