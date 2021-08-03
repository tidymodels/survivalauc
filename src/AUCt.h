#ifndef _AUCT_H_
#define _AUCT_H_

void AUCt(int *n, double *Y, double *D, double *Z, int *k, double *threshold, int *m, double *uniqY,
          double **sens, double **spec, double *auc);
void coxph_surv(int *n, double *Y, double *D, double *Z, int *m, double *uniqY, int *locY,
                double **SurvEst, double *ESurvEst);
void ESurv(int n, double *Z, int m, double *uniqY, double *Lambda, double beta, double **SurvEst, double *ESurvEst);
void sens_spec(int n, int m, int k, double *Z, double *threshold, double **SurvEst, double *ESurvEst,
               double **sens, double **spec);
void AUC(int n, int m, double *Z, double **SurvEst, double *ESurvEst, double *auc);

#endif
