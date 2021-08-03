/*
 *  CoxPH.c
 *
 *  Fuction: Implementing the nonparametric maximum likelihood estimation (NPMLE) of the baseline
 *           cumulative hazard function and the regression parameters under the Cox model
 *           for right-censored data. The NPMLE of the regression coefficients is the same as
 *           the maximum partial likelihood estimators and the NPMLE of the baseline cumulative
 *           hazard function is the Brewlow estimator.
 *
 *  Syntax:
 *      void coxph(int n, double *Y, double *Delta, int p, double **X, double *beta,
 *               double *lambda, double *Lambda, int *L, double *uniqY, int *locY,
 *               double *betase, double **cov)
 *  Input arguments:
 *      n = sample size
 *      Y = survival/censoring time, n by 1 vector
 *      Delta = censoring indicator taking value 1 if survival time is observed
 *          and 0 otherwise, n by 1 vector
 *      p = number of covariates
 *      X = n by p design matrix
 *
 *  Output arguments:
 *      L = number of distinct observed failure times
 *      uniqY = distinct observed failure times in ascending order; m by 1 vector
 *      locY = n by 1 vector; locY[i] = j, where j is the largest integer such that uniqY[j] <= Y[i]
 *      beta = NPMLE of regression coefficients; p by 1 vector
 *      lambda = baseline hazard estimate; m by 1 vector; ith element corresponds to the
 *          jump size of baseline cumulative hazard function at uniqY[i]
 *      Lambda = baseline cumulative hazard estimate; m by 1 vector; ith element corresponds
 *          to the baseline cumulative hazard estimate at uniqY[i]
 *      betase = standard error estimate of betahat; p by 1 vector
 *      cov = (L+p) by (L+p) matrix; inverse of the negative second derivatives of
 *          the log-likelihood (treating jump sizes of baseline cumulative hazard at the
 *          observed failure times and the regression coefficients as unknown parameters)
 *
 *  Author: Guoqing Diao
 *  Last update date: 06/25/2021
 *
 */

#include "CoxPH.h"
#include "nrutil.h"
#include <math.h>

/*
 *  Subroutine: locate_pos()
 *
 *  Fuction: Given a sequence of numbers in ascending order a[start:end], locate the interval
 *          that covers x, where x >= a[start]. Returns to an integer start <= l <= end such that
 *          a[l] <= x < a[l+1]. Used in pre-processing right-censored survival data
 *
 *  Input arguments:
 *      a = a vector
 *      start = starting location
 *      end = ending location
 *      x = a double number
 *
 *  Return value: an integer start <= l <= end such that a[l] <= x < a[l+1].
 *
 */

int locate_pos(double *a, int start, int end, double x)
{
        int i,bn,loc,mina;

        bn = (int)((end+start)/2);

        if (x >= a[end])
        {
                return end;
        }
        else if (x >= a[bn] && x < a[bn+1])
        {
                return(bn);
        }
        else if (x < a[bn])
        {
                return locate_pos(a, start, bn-1,x);
        }
        else if (x >= a[bn+1])
        {
                return locate_pos(a, bn+1, end, x);
        }

    return 0;
}

/*
 *  Subroutine: DataM()
 *
 *  Fuction: Pre-processing right-censored survival data
 *
 *  Input arguments:
 *      n = sample size
 *      Y = survival/censoring time, n by 1 vector
 *      Delta = censoring indicator taking value 1 if survival time is observed
 *          and 0 otherwise, n by 1 vector
 *
 *  Output arguments:
 *      L = number of distinct observed failure times
 *      uniqY = distinct observed failure times in ascending order; m by 1 vector
 *      locY = n by 1 vector; locY[i] = j, where j is the largest integer such that uniqY[j] <= Y[i]
 *
 *
 */
void DataM(int n, double *Y, double *D, int *L, double *uniqY, int *locY)
{
        int i,j,k,l,*d;
        double *ft,*delta;

        ft = dvector(1,n); d = ivector(1,n);
        delta = dvector(1,n);

        for (i=1; i<=n; i++)
        {
                ft[i] = Y[i];
                delta[i] = D[i];
        }

        sort2(n,ft,delta);
        i = 0; *L = 0;

        for (j=1; j<=n; j++)
                d[j] = 0;
        for (j=1; j<=n; j++)
        {
                i++;
                if (delta[j]==1) d[i]++;
                while (j<n && ft[j+1]==ft[j])
                {
                        if (delta[j+1]==1)
                                d[i]++;
                        j++;
                }
                if (d[i]==0)
                        i--;
                else
                {
                        uniqY[i] = ft[j];
                }
        }
       *L = i;

        for (i=1; i<=n; i++)
        {
                if (Y[i] < uniqY[1] && D[i] == 0) locY[i] = 0;
                else if (Y[i] < uniqY[1] && D[i] == -1) locY[i] = 1;
                else if (Y[i] >= uniqY[*L] && D[i] == -1) locY[i] = 0;
                else if (D[i] == -1) locY[i]= locate_pos(uniqY,1,*L,Y[i]) + 1;
                else locY[i] = locate_pos(uniqY,1,*L,Y[i]);
        }

        free_dvector(ft,1,n); free_dvector(delta,1,n);
        free_ivector(d,1,n);
}

// Define external variables
int coxph_En, coxph_Ep, *coxph_ElocY, coxph_EL;
double *coxph_EDelta, *coxph_EY, *coxph_EuniqY, **coxph_EX, *coxph_Ebeta,
	*coxph_Elambda, *coxph_ELambda, *coxph_Ed1beta,
	*coxph_Ed1lambda, *coxph_theta, *coxph_d1theta,
	**coxph_d2theta, **coxph_cov;

// Allocate memory to external variables
void coxph_alloc_memory(long n, long p)
{
        coxph_ElocY = ivector(1, n);
        coxph_EDelta = dvector(1, n);
        coxph_EY = dvector(1, n);
        coxph_EuniqY = dvector(1, n);
        coxph_EX = dmatrix(1, n, 1, p);
        coxph_Ebeta = dvector(1, p);
        coxph_Elambda = dvector(1, n);
        coxph_ELambda = dvector(1, n);

        coxph_Ed1beta = dvector(1, p);
        coxph_Ed1lambda = dvector(1, n);
        coxph_theta = dvector(1, p+n);
        coxph_d1theta = dvector(1, p+n);
        coxph_d2theta = dmatrix(1, p+n, 1, p+n);
        coxph_cov = dmatrix(1, p+n, 1, p+n);
}
// Free memory for external variables
void coxph_free_memory(long n, long p)
{
        free_ivector(coxph_ElocY, 1, n);
        free_dvector(coxph_EDelta, 1, n);
        free_dvector(coxph_EY, 1, n);
        free_dvector(coxph_EuniqY, 1, n);
        free_dmatrix(coxph_EX, 1, n, 1, p);
        free_dvector(coxph_Ebeta, 1, p);
        free_dvector(coxph_Elambda, 1, n);
        free_dvector(coxph_ELambda, 1, n);

        free_dvector(coxph_Ed1beta, 1, p);
        free_dvector(coxph_Ed1lambda, 1, n);
        free_dvector(coxph_theta, 1, p+n);
        free_dvector(coxph_d1theta, 1, p+n);
        free_dmatrix(coxph_d2theta, 1, p+n, 1, p+n);
        free_dmatrix(coxph_cov, 1, p+n, 1, p+n);

}

/*
 *  Subroutine: coxph_logL()
 *
 *  Fuction: calculate the (nonparametric) log-likelihood function of unknown parameters,
 *              including regression coefficients and jump sizes of the baseline cumulative
 *              hazard at the observed failure times
 *
 *  Input arguments:
 *      n = sample size
 *      Delta = censoring indicator taking value 1 if survival time is observed
 *          and 0 otherwise, n by 1 vector
 *      p = number of covariates
 *      X = n by p design matrix
 *      L = number of distinct observed failure times
 *      locY = n by 1 vector; locY[i] = j, where j is the largest integer such that uniqY[j] <= Y[i]
 *      beta = NPMLE of regression coefficients; p by 1 vector
 *      lambda = baseline hazard estimate; m by 1 vector; ith element corresponds to the
 *          jump size of baseline cumulative hazard function at uniqY[i]
 *
 *  Return value: log-likelihood
 *
 */
double coxph_logL(int n, int p, int *locY, double *Delta, double **X, double *beta,
        int L, double *lambda)
{
        int i, j, k, l;
        double *Lambda, betaX, expbetaX, sum;

        Lambda = dvector(1, L);

        Lambda[1] = lambda[1];
        for (j=2; j<=L; j++)
                Lambda[j] = Lambda[j-1] + lambda[j];

        sum = 0;
        for (i=1; i<=n; i++)
        {
                betaX = 0;
                for (j=1; j<=p; j++)
                        betaX += beta[j]*X[i][j];
                expbetaX = exp(betaX);

		if (locY[i] > 0)
		sum += Delta[i]*(betaX + log(lambda[locY[i]])) - Lambda[locY[i]]*expbetaX;
        }

        free_dvector(Lambda, 1, L);
        return sum;
}

/*
 *  Subroutine: coxph_d1logL()
 *
 *  Fuction: calculate the first derivatives of the (nonparametric) log-likelihood function
 *              with respect to the regression coefficients and the jump sizes of the baseline
 *              cumulative hazard at the observed failure times
 *
 *  Input arguments:
 *      n = sample size
 *      Delta = censoring indicator taking value 1 if survival time is observed
 *          and 0 otherwise, n by 1 vector
 *      p = number of covariates
 *      X = n by p design matrix
 *      L = number of distinct observed failure times
 *      locY = n by 1 vector; locY[i] = j, where j is the largest integer such that uniqY[j] <= Y[i]
 *      beta = NPMLE of regression coefficients; p by 1 vector
 *      lambda = baseline hazard estimate; m by 1 vector; ith element corresponds to the
 *          jump size of baseline cumulative hazard function at uniqY[i]
 *
 *  Output arguments:
 *      d1beta = first derivatives with respect to beta; p by 1 vector
 *      d1lambda = first derivatives with respect to the jump sizes of the baseline
 *          cumulative hazard at the observed failure times; L by 1 vector
 *
 */
void coxph_d1logL(int n, int p, int *locY, double *Delta, double **X, double *beta,
        int L, double *lambda, double *d1beta, double *d1lambda)
{
        int i, j, k, l;
        double *Lambda, betaX, expbetaX, tmp;

        Lambda = dvector(1, L);

        Lambda[1] = lambda[1];
        for (j=2; j<=L; j++)
                Lambda[j] = Lambda[j-1] + lambda[j];

        for (j=1; j<=p; j++)
                d1beta[j] = 0;
        for (j=1; j<=L; j++)
                d1lambda[j] = 0;
        for (i=1; i<=n; i++)
        {
		if (locY[i] > 0) {
                betaX = 0;
                for (j=1; j<=p; j++)
                        betaX += beta[j]*X[i][j];
                expbetaX = exp(betaX);

                for (j=1; j<=p; j++)
                        d1beta[j] += (Delta[i]-Lambda[locY[i]]*expbetaX)*X[i][j];

                for (j=1; j<=locY[i]; j++)
                        d1lambda[j] -= expbetaX;
                if (Delta[i] == 1)
                        d1lambda[locY[i]] += 1/lambda[locY[i]];
		}
        }
        free_dvector(Lambda, 1, L);
}

/*
 *  Subroutine: coxph_d2logL()
 *
 *  Fuction: calculate the negative second derivatives of the (nonparametric)
 *              log-likelihood function with respect to the regression coefficients
 *              and the jump sizes of the baseline cumulative hazard at the observed
 *              failure times; the inverse is the estimate of the covariance matrix
 *              of the NPMLEs
 *
 *  Input arguments:
 *      n = sample size
 *      Delta = censoring indicator taking value 1 if survival time is observed
 *          and 0 otherwise, n by 1 vector
 *      p = number of covariates
 *      X = n by p design matrix
 *      L = number of distinct observed failure times
 *      locY = n by 1 vector; locY[i] = j, where j is the largest integer such that uniqY[j] <= Y[i]
 *      beta = NPMLE of regression coefficients; p by 1 vector
 *      lambda = baseline hazard estimate; m by 1 vector; ith element corresponds to the
 *          jump size of baseline cumulative hazard function at uniqY[i]
 *
 *  Output arguments:
 *      d2theta = negative second derivatives with respect to beta and the jump sizes of
 *          the baseline cumulative hazard at the observed failure times; (p+L) by (p+L) matrix
 *
 */
void coxph_d2logL(int n, int p, int *locY, double *Delta, double **X, double *beta,
        int L, double *lambda, double **d2theta)
{
        int i, j, k, l;
        double *Lambda, betaX, expbetaX, tmp;

        Lambda = dvector(1, L);

        Lambda[1] = lambda[1];
        for (j=2; j<=L; j++)
                Lambda[j] = Lambda[j-1] + lambda[j];

        for (j=1; j<=p+L; j++)
                for (k=1; k<=p+L; k++)
                        d2theta[j][k] = 0;

        for (i=1; i<=n; i++)
        {
		if (locY[i] > 0)
		{
                betaX = 0;
                for (j=1; j<=p; j++)
                        betaX += beta[j]*X[i][j];
                expbetaX = exp(betaX);

                // d beta d beta
                tmp = Lambda[locY[i]]*expbetaX;
                for (j=1; j<=p; j++)
                        for (k=1; k<=j; k++)
                                d2theta[j][k] += tmp*X[i][j]*X[i][k];
                // d beta d lambda
                for (j=1; j<=p; j++)
                        for (l=1; l<=locY[i]; l++)
                                d2theta[p+l][j] += expbetaX*X[i][j];

                // d lambda d lambda
                d2theta[p+locY[i]][p+locY[i]] += Delta[i]/
                        pow(lambda[locY[i]],2.0);
		}
        }

        for (j=1; j<=p+L; j++)
                for (k=1; k<=j; k++)
        		d2theta[k][j] = d2theta[j][k];

        free_dvector(Lambda, 1, L);
}

// function used in the quasi-Newton algorithm
// negative log-likelihood
double coxph_ElogL(double *theta)
{
        int i, j, k;
        double sum;

        for (j=1; j<=coxph_Ep; j++)
                coxph_Ebeta[j] = theta[j];
        for (j=1; j<=coxph_EL; j++)
                coxph_Elambda[j] = exp(theta[coxph_Ep+j]);

        sum = -coxph_logL(coxph_En, coxph_Ep, coxph_ElocY, coxph_EDelta, coxph_EX,
                coxph_Ebeta, coxph_EL, coxph_Elambda);
        return sum;
}

// function used in the quasi-Newton algorithm
// first derivatives of the negative log-likelihood function
void coxph_Ed1logL(double *theta, double *d1theta)
{
        int i, j, k;

        for (j=1; j<=coxph_Ep; j++)
                coxph_Ebeta[j] = theta[j];
        for (j=1; j<=coxph_EL; j++)
                coxph_Elambda[j] = exp(theta[coxph_Ep+j]);

        coxph_d1logL(coxph_En, coxph_Ep, coxph_ElocY, coxph_EDelta, coxph_EX,
                coxph_Ebeta, coxph_EL, coxph_Elambda,
                coxph_Ed1beta, coxph_Ed1lambda);
        for (j=1; j<=coxph_Ep; j++)
        {
                d1theta[j] = -coxph_Ed1beta[j];
        }
        for (j=1; j<=coxph_EL; j++)
                d1theta[coxph_Ep+j] = -coxph_Ed1lambda[j]*coxph_Elambda[j];
}

// second derivatives of the negative log-likelihood function
void coxph_Ed2logL(double *theta, double **d2theta)
{
        int i, j, k;

        for (j=1; j<=coxph_Ep; j++)
                coxph_Ebeta[j] = theta[j];
        for (j=1; j<=coxph_EL; j++)
                coxph_Elambda[j] = exp(theta[coxph_Ep+j]);

        coxph_d2logL(coxph_En, coxph_Ep, coxph_ElocY, coxph_EDelta, coxph_EX, coxph_Ebeta,
                coxph_EL, coxph_Elambda, d2theta);
}


void coxph(int n, double *Y, double *Delta, int p, double **X, double *beta, double *lambda,
	double *Lambda, int *L, double *uniqY, int *locY, double *betase, double **cov)
{
	int i, j, k, iter;
	double f1;

	// allocate memory
	coxph_alloc_memory(n, p);

	// Data
	coxph_En = n;
	coxph_Ep = p;
	for (i=1; i<=n; i++)
	{
		coxph_EY[i] = Y[i];
		coxph_EDelta[i] = Delta[i];
		for (j=1; j<=p; j++)
			coxph_EX[i][j] = X[i][j];
	}
	// Data Management
    DataM(coxph_En, coxph_EY, coxph_EDelta, &coxph_EL, coxph_EuniqY, coxph_ElocY);

	// initial values for theta
	for (j=1; j<=coxph_Ep; j++)
		coxph_theta[j] = 0;
	for (j=1; j<=coxph_EL; j++)
		coxph_theta[coxph_Ep+j] = log(1.0/coxph_EL);

    // call the quasi-Newton algorithm
	dfpmin(coxph_theta, coxph_Ep+coxph_EL, 1.0e-12, &iter,
		&f1, coxph_ElogL, coxph_Ed1logL);


	// beta
	for (j=1; j<=p; j++)
		beta[j] = coxph_theta[j];
	// lambda, Lambda
	lambda[1] = exp(coxph_theta[p+1]); Lambda[1] = lambda[1];
	for (j=2; j<=coxph_EL; j++)
	{
		lambda[j] = exp(coxph_theta[p+j]);
		Lambda[j] = Lambda[j-1] + lambda[j];
	}
	// L
	*L = coxph_EL;
	// uniqY
	for (j=1; j<=coxph_EL; j++)
		uniqY[j] = coxph_EuniqY[j];
	// locY
	for (i=1; i<=n; i++)
		locY[i] = coxph_ElocY[i];
	// d2theta
	coxph_Ed2logL(coxph_theta, coxph_d2theta);
	// cov
	dmatrix_inv(coxph_d2theta, cov, coxph_Ep+coxph_EL);
	// betase
	for (j=1; j<=p; j++)
		betase[j] = sqrt(cov[j][j]);

	coxph_free_memory(n, p);
}

// Given the NPMLEs of beta and jump sizes of baseline cumulative hazard, calculate
// the NPMLE of the baseline cumulative hazard and the standard error estimates
void coxph_Lambdat(int p, int L, double *uniqY, double *Lambda,
        double **cov, int nt, double *t, double *Lambdat,
        double *Lambdat_se)
{
        int i, j, k, l;
        double tmp;

        for (l=1; l<=nt; l++)
        {
                k = locate_pos(uniqY, 1, L, t[l]);
                tmp = 0;
                for (i=1; i<=k; i++)
                        for (j=1; j<=k; j++)
                                tmp += cov[p+i][p+j];

                Lambdat[l] = Lambda[k];
                Lambdat_se[l] = sqrt(tmp);
        }
}

