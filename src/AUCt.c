/*
 *  AUCt.c
 *
 *  Fuction: Implementing the alternative estimators of the time-dependent area
 *           under the ROC curve, time-dependent sensitivity and specificity,
 *           and time-dependent ROC curves. The corresponding estimators are
 *           presented in equations (9)-(11) in Chambless and Diao (Statistics
 *           in Medicine, 2006; 25: 3474-3486).
 *
 *  Input arguments:
 *      n = sample size
 *      Y = survival/censoring time, n by 1 vector
 *      D = censoring indicator taking value 1 if survival time is observed
 *          and 0 otherwise, n by 1 vector
 *      Z = risk score obtained from a risk prediction model, n by 1 vector
 *      k = number of thresholds used in estimating sensitivity and specificity
 *      threshold = k by 1 vector; thresholds used in estimating sensitivity and specificity
 *
 *  Output arguments:
 *      m = number of distinct observed failure times
 *      uniqY = distinct observed failure times in ascending order; m by 1 vector
 *      sens = m by k matrix, the (i,j)th element corresponds to the estimated sensitivity at time point
 *          uniqY[i] with threshold thresold[j]
 *      spec = m by k matrix, the (i,j)th element corresponds to the estimated specificity at time point
 *          uniqY[i] with threshold threshold[j]
 *      auc = m by 1 vector, the ith element corresponds to the estimated AUC at time point
 *          uniqY[i]
 *
 *
 *  Author: Guoqing Diao
 *  Last update date: 06/25/2021
 *
 */

#include "AUCt.h"
#include "CoxPH.h"
#include "nrutil.h"
#include <math.h>

void AUCt(int *n, double *Y, double *D, double *Z, int *k, double *threshold, int *m, double *uniqY,
          double **sens, double **spec, double *auc)
{
    int p, *locY;
    double **SurvEst, *ESurvEst;

    p = 1; // number of covariates
    locY = ivector(1, *n);

    // Data processing before fitting survival models
    // Output:
    //  m = number of distinct observed failure times
    //  uniqY = distinct observed failure times in an ascending order, m by 1 vector
    //  locY = n by 1 vector; locY[i] = j, where j is the largest integer such that uniqY[j] <= Y[i]
    DataM(*n, Y, D, m, uniqY, locY);

    // Fit a Cox model with risk score Z as covariates to obtain the survival function estimates
    // Output:
    //  SurvEst = m by n matrix, the [i,j] element corresponds to the survival function estimate at
    //              the ith ordered distinct observed failure time (uniqY[i]) given the jth subject's
    //              risk score (Z[j]); SurvEst[i,j] = exp(-Lambda[i]*exp(beta*Z[j]))
    //  ESurvEst = m by 1 vector, empirical estimate of E(S(t|Z)), the ith element is the sample average
    //              of the jth row of SurvEst (SurvEst[i,])
    // Remark: other survival models may also be used
    SurvEst = dmatrix(1, *m, 1, *n);
    ESurvEst = dvector(1, *m);
    coxph_surv(n, Y, D, Z, m, uniqY, locY, SurvEst, ESurvEst);

    // estimate sensitivity and specificity at each distinct observed failure time with a set of thresholds
    sens_spec(*n, *m, *k, Z, threshold, SurvEst, ESurvEst, sens, spec);

    // estimate AUC at each distinct observed failure time
    AUC(*n, *m, Z, SurvEst, ESurvEst, auc);
}

/*
 *  Subroutine: coxph_surv()
 *
 *  Fuction: Obtain the conditional survival fucntion estimates S(t|Z) based on the Cox model and
 *              the empirical estimate of E(S(t|Z))
 *
 *  Input arguments:
 *      n = sample size
 *      Y = survival/censoring time, n by 1 vector
 *      D = censoring indicator taking value 1 if survival time is observed
 *          and 0 otherwise, n by 1 vector
 *      Z = risk score obtained from a risk prediction model, n by 1 vector
 *      m = number of distinct observed failure times
 *      uniqY = distinct observed failure times in ascending order, m by 1 vector
 *      locY = n by 1 vector; locY[i] = j, where j is the largest integer such that uniqY[j] <= Y[i]
 *
 *  Output arguments:
 *      SurvEst = m by n matrix, the [i,j] element corresponds to the survival function estimate at
 *              the ith ordered distinct observed failure time (uniqY[i]) given the jth subject's
 *             risk score (Z[j]); SurvEst[i,j] = exp(-Lambda[i]*exp(beta*Z[j]))
 *      ESurvEst = m by 1 vector, empirical estimate of E(S(t|Z)), the ith element is the sample average
 *              of the ith row of SurvEst (SurvEst[i,])
 *
 *  Remark: other survival models can also be use to estimate the conditional survival function S(t|Z)
 */
void coxph_surv(int *n, double *Y, double *D, double *Z, int *m, double *uniqY, int *locY,
                double **SurvEst, double *ESurvEst)
{
    int i, p;
    double **tmpZ, *beta, *betase, **cov, *lambda, *Lambda;

    p = 1; // number of covariates
    tmpZ = dmatrix(1, *n, 1, p); // convert Z to a n by 1 matrix
    beta = dvector(1, p); // regression coefficient in the Cox model
    for (i=1; i<=*n; i++)
        tmpZ[i][1] = Z[i];
    betase = dvector(1, p); // standard error estimate of betahat
    cov = dmatrix(1, *n+p, 1, *n+p); // covariance matrix estimate
    lambda = dvector(1, *n);
    Lambda = dvector(1, *n);


    // Fit a Cox model with risk score Z as the covariate to obtain the survival function estimates
    // Output:
    //  lambda = baseline hazard estimate; m by 1 vector; ith element corresponds to the jump size of
    //      baseline cumulative hazard function at uniqY[i]
    //  Lambda = baseline cumulative hazard estimate; m by 1 vector; ith element corresponds to the baseline
    //      cumulative hazard estimate at uniqY[i]
    //  beta = regression coefficient estimate
    //  betase = standard error estimate of betahat
    //  cov = (m+p) by (m+p) matrix, inverse of the negative second derivatives of the log-likelihood
    coxph(*n, Y, D, p, tmpZ, beta, lambda, Lambda, m, uniqY, locY, betase, cov);

    // call the ESurv subroutine to calculate the survival estimates at each distinct observed failure time
    // for each subject and the empirical estimate of E(S(t|Z))
    ESurv(*n, Z, *m, uniqY, Lambda, beta[1], SurvEst, ESurvEst);
}

/*
 *  Subroutine: ESurv()
 *
 *  Fuction: Obtain the empirical estimate of E(S(t|Z)), where S(t|Z) is the survival function at
 *           time point t given risk score Z.
 *
 *  Input arguments:
 *      n = sample size
 *      Z = risk score; n by 1 vector
 *      m = number of distinct observed failure times
 *      uniqY = distinct observed failure times in ascending order, m by 1 vector
 *      Lambda = baseline cumulative hazard function estimate (from the Cox model), m by 1 vector
 *      beta = regression coefficient estimate of the effect of risk score Z (from the Cox model)
 *
 *  Output arguments:
 *      SurvEst = m by n matrix, the [i,j] element corresponds to the survival function estimate at
 *              the ith ordered distinct observed failure time (uniqY[i]) given the jth subject's
 *             risk score (Z[j]); SurvEst[i,j] = exp(-Lambda[i]*exp(beta*Z[j]))
 *      ESurvEst = m by 1 vector, empirical estimate of E(S(t|Z)), the ith element is the sample average
 *              of the jth row of SurvEst (SurvEst[i,])
 *
 */
void ESurv(int n, double *Z, int m, double *uniqY, double *Lambda, double beta, double **SurvEst, double *ESurvEst)
{
    int i, j;

    for (i=1; i<=m; i++)
    {
        ESurvEst[i] = 0;
        for (j=1; j<=n; j++)
        {
            SurvEst[i][j] = exp(-Lambda[i]*exp(beta*Z[j]));
            ESurvEst[i] += SurvEst[i][j];
        }
        ESurvEst[i] /= n;
    }
}

/*
 *  Subroutine: sens_spec()
 *
 *  Fuction: Given the survival estimates at each distinct observed failure time for each subject,
 *              estimate the sensitivity and specificity at each distinct observed failure time for
 *              a vector of thresholds.
 *
 *  Input arguments:
 *      n = sample size
 *      m = number of distinct observed failure times
 *      k = number of thresholds
 *      Z = risk score obtained from a risk prediction model, n by 1 vector
 *      threshold = k by 1 vector
 *      SurvEst = m by n matrix, the [i,j] element corresponds to the survival function estimate at
 *              the ith ordered distinct observed failure time (uniqY[i]) given the jth subject's
 *             risk score (Z[j]); SurvEst[i,j] = exp(-Lambda[i]*exp(beta*Z[j]))
 *      ESurvEst = m by 1 vector, empirical estimate of E(S(t|Z)), the ith element is the sample average
 *              of the jth row of SurvEst (SurvEst[i,])
 *
 *  Output arguments:
 *      sens = m by k matrix, the (i,j)th element corresponds to the estimated sensitivity at time point
 *          uniqY[i] with threshold thresold[j]
 *      spec = m by k matrix, the (i,j)th element corresponds to the estimated specificity at time point
 *          uniqY[i] with threshold threshold[j]
 *
 */
void sens_spec(int n, int m, int k, double *Z, double *threshold, double **SurvEst, double *ESurvEst,
               double **sens, double **spec)
{
    int i, j, l;
    double sens_num, spec_num;

    for (i=1; i<=m; i++)
        for (j=1; j<=k; j++)
        {
            // time point: uniqY[i]
            // threshold: threshold[j]
            // calculate the numerator terms in equations (10) and (11) of Chambless and Diao (2006)
            sens_num = 0; spec_num = 0;
            for (l=1; l<=n; l++)
            {
                if (Z[l] > threshold[j])
                    sens_num += 1 - SurvEst[i][l];
                else spec_num += SurvEst[i][l];
            }
            sens_num /= n;
            spec_num /= n;
            sens[i][j] = sens_num/(1-ESurvEst[i]);
            spec[i][j] = spec_num/ESurvEst[i];
        }
}

/*
 *  Subroutine: AUC()
 *
 *  Fuction: Given the survival estimates at each distinct observed failure time for each subject,
 *              estimate the AUC at each distinct observed failure time.
 *
 *  Input arguments:
 *      n = sample size
 *      m = number of distinct observed failure times
 *      Z = risk score obtained from a risk prediction model, n by 1 vector
 *      SurvEst = m by n matrix, the [i,j] element corresponds to the survival function estimate at
 *              the ith ordered distinct observed failure time (uniqY[i]) given the jth subject's
 *             risk score (Z[j]); SurvEst[i,j] = exp(-Lambda[i]*exp(beta*Z[j]))
 *      ESurvEst = m by 1 vector, empirical estimate of E(S(t|Z)), the ith element is the sample average
 *              of the jth row of SurvEst (SurvEst[i,])
 *
 *  Output arguments:
 *      auc = m by 1 vector, the ith element corresponds to the estimated AUC at time point
 *          uniqY[i]
 *
 */
void AUC(int n, int m, double *Z, double **SurvEst, double *ESurvEst, double *auc)
{
    int i, j, l;
    double auc_num;

    for (i=1; i<=m; i++)
    {
        // time point: uniqY[i]
        // numerator term in equation (9) in Chambless and Diao (2006)
        auc_num = 0;
        for (j=1; j<=n-1; j++)
            for (l=j+1; l<=n; l++)
        {
            if (Z[j] > Z[l]) auc_num += (1 - SurvEst[i][j])*SurvEst[i][l];
            else if (Z[j] < Z[l]) auc_num += SurvEst[i][j]*(1 - SurvEst[i][l]);
        }
        auc_num /= n*(n-1);

        auc[i] = auc_num/(ESurvEst[i]*(1-ESurvEst[i]));
    }
}
