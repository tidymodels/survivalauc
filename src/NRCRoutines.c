/*
 Routines from Numerical Recipes in C
 */
#define NR_END 1
#define FREE_ARG char*
int ludcmp_flag, flag;

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
        fprintf(stderr,"Numerical Recipes run-time error...\n");
        fprintf(stderr,"%s\n",error_text);
        fprintf(stderr,"...now exiting to system...\n");
        exit(1);
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
        float *v;

        v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
        if (!v) nrerror("allocation failure in vector()");
        return v-nl+NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;

        v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
        if (!v) nrerror("allocation failure in ivector()");
        return v-nl+NR_END;
}

unsigned char *cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
        unsigned char *v;

        v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
        if (!v) nrerror("allocation failure in cvector()");
        return v-nl+NR_END;
}
unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
        unsigned long *v;

        v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
        if (!v) nrerror("allocation failure in lvector()");
        return v-nl+NR_END;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
        double *v;

        v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
        if (!v) nrerror("allocation failure in dvector()");
        return v-nl+NR_END;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        double **m;

        /* allocate pointers to rows */
        m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
        if (!m) nrerror("allocation failure 1 in matrix()");
        m += NR_END;
        m -= nrl;

        /* allocate rows and set pointers to them */
        m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
        if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
        m[nrl] += NR_END;
        m[nrl] -= ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}
int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        int **m;

        /* allocate pointers to rows */
        m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
        if (!m) nrerror("allocation failure 1 in matrix()");
        m += NR_END;
        m -= nrl;


        /* allocate rows and set pointers to them */
        m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
        if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
        m[nrl] += NR_END;
        m[nrl] -= ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}

float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
        long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
        long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
        float **m;

        /* allocate array of pointers to rows */
        m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
        if (!m) nrerror("allocation failure in submatrix()");
        m += NR_END;
        m -= newrl;

        /* set pointers to rows */
        for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}

double **subdmatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch,
        long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
        long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
        double **m;

        /* allocate array of pointers to rows */
        m=(double **) malloc((size_t) ((nrow+NR_END)*sizeof(double*)));
        if (!m) nrerror("allocation failure in submatrix()");
        m += NR_END;
        m -= newrl;

        /* set pointers to rows */
        for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}

float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
        long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
        float **m;

        /* allocate pointers to rows */
        m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
        if (!m) nrerror("allocation failure in convert_matrix()");
        m += NR_END;
        m -= nrl;

        /* set pointers to rows */
        m[nrl]=a-ncl;
        for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
        /* return pointer to array of pointers to rows */
        return m;
}
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
        long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
        float ***t;

        /* allocate pointers to pointers to rows */
        t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
        if (!t) nrerror("allocation failure 1 in f3tensor()");
        t += NR_END;
        t -= nrl;

        /* allocate pointers to rows and set pointers to them */
        t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
        if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
        t[nrl] += NR_END;
        t[nrl] -= ncl;

        /* allocate rows and set pointers to them */
        t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
        if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
        t[nrl][ncl] += NR_END;
        t[nrl][ncl] -= ndl;

        for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
        for(i=nrl+1;i<=nrh;i++) {
                t[i]=t[i-1]+ncol;
                t[i][ncl]=t[i-1][ncl]+ncol*ndep;
                for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
        }

        /* return pointer to array of pointers to rows */
        return t;
}
void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(unsigned char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}

void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
        free((FREE_ARG) (b+nrl-NR_END));
}

void free_subdmatrix(double **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
        free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a matrix allocated by convert_matrix() */
{
        free((FREE_ARG) (b+nrl-NR_END));
}

void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
        long ndl, long ndh)
/* free a float f3tensor allocated by f3tensor() */
{
        free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
        free((FREE_ARG) (t[nrl]+ncl-NR_END));
        free((FREE_ARG) (t+nrl-NR_END));
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5.){2p491&.#@#Q0JS[. */

void lubksb(double **a, int n, int *indx, double b[])
{
        int i,ii=0,ip,j;
        double sum;

        for (i=1;i<=n;i++) {
                ip=indx[i];
                sum=b[ip];
                b[ip]=b[i];
                if (ii)
                        for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
                else if (sum) ii=i;
                b[i]=sum;
        }
        for (i=n;i>=1;i--) {
                sum=b[i];
                for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
                b[i]=sum/a[i][i];
        }
}
void ludcmp(double **a, int n, int *indx, double *d)
{
        int i,imax,j,k;
        double big,dum,sum,temp;
        double *vv, TINY = 1.0e-20;

        vv=dvector(1,n);
        *d=1.0;
        for (i=1;i<=n;i++) {
                big=0.0;
                for (j=1;j<=n;j++)
                        if ((temp=fabs(a[i][j])) > big) big=temp;
                if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
                vv[i]=1.0/big;
        }
        for (j=1;j<=n;j++) {
                for (i=1;i<j;i++) {
                        sum=a[i][j];
                        for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
                        a[i][j]=sum;
                }
                big=0.0;
                for (i=j;i<=n;i++) {
                        sum=a[i][j];
                        for (k=1;k<j;k++)
                                sum -= a[i][k]*a[k][j];
                        a[i][j]=sum;
                        if ( (dum=vv[i]*fabs(sum)) >= big) {
                                big=dum;
                                imax=i;
                        }
                }
                if (j != imax) {
                        for (k=1;k<=n;k++) {
                                dum=a[imax][k];
                                a[imax][k]=a[j][k];
                                a[j][k]=dum;
                        }
                        *d = -(*d);
                        vv[imax]=vv[j];
                }
                indx[j]=imax;
                if (a[j][j] == 0.0) a[j][j]=TINY;
                if (j != n) {
                        dum=1.0/(a[j][j]);
                        for (i=j+1;i<=n;i++) a[i][j] *= dum;
                }
        }
        free_dvector(vv,1,n);
}

#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349

double ran1(idum)
long *idum;
{
        static long ix1,ix2,ix3;
        static double r[98];
        double temp;
        static int iff=0;
        int j;
        void nrerror();

        if (*idum < 0 || iff == 0) {
                iff=1;
                ix1=(IC1-(*idum)) % M1;
                ix1=(IA1*ix1+IC1) % M1;
                ix2=ix1 % M2;
                ix1=(IA1*ix1+IC1) % M1;
                ix3=ix1 % M3;
                for (j=1;j<=97;j++) {
                        ix1=(IA1*ix1+IC1) % M1;
                        ix2=(IA2*ix2+IC2) % M2;
                        r[j]=(ix1+ix2*RM2)*RM1;
                }
                *idum=1;
        }
        ix1=(IA1*ix1+IC1) % M1;
        ix2=(IA2*ix2+IC2) % M2;
        ix3=(IA3*ix3+IC3) % M3;
        j=1 + ((97*ix3)/M3);
        if (j > 97 || j < 1) nrerror("RAN1: This cannot happen.");
        temp=r[j];
        r[j]=(ix1+ix2*RM2)*RM1;
        return temp;
}

#undef M1
#undef IA1
#undef IC1
#undef RM1
#undef M2
#undef IA2
#undef IC2
#undef RM2
#undef M3
#undef IA3
#undef IC3

double gasdev(idum)
long *idum;
{
        static int iset=0;
        static double gset;
        double fac,r,v1,v2;
        double ran2();

        if  (iset == 0) {
                do {
                        v1=2.0*ran2(idum)-1.0;
                        v2=2.0*ran2(idum)-1.0;
                        r=v1*v1+v2*v2;
                } while (r >= 1.0);
                fac=sqrt(-2.0*log(r)/r);
                gset=v1*fac;
                iset=1;
                return v2*fac;
        } else {
                iset=0;
                return gset;
        }
}

double gamdev(int ia, long *idum)
{
        double ran2(long *idum);
        void nrerror(char error_text[]);
        int j;
        double am,e,s,v1,v2,x,y;

        if (ia < 1) nrerror("Error in routine gamdev");
        if (ia < 6) {
                x=1.0;
                for (j=1;j<=ia;j++) x *= ran2(idum);
                x = -log(x);
        } else {
                do {
                        do {
                                do {
                                        v1=2.0*ran2(idum)-1.0;
                                        v2=2.0*ran2(idum)-1.0;
                                } while (v1*v1+v2*v2 > 1.0);
                                y=v2/v1;
                                am=ia-1;
                                s=sqrt(2.0*am+1.0);
                                x=s*y+am;
                        } while (x <= 0.0);
                        e=(1.0+y*y)*exp(am*log(x/am)-s*y);
                } while (ran2(idum) > e);
        }
        return x;
}

void tdev(int df, int dim, double *trannum, long *idum)
{
	int i,j;
	double Y;

	Y = gamdev(df/2.0,idum)*2; Y = sqrt(Y/df);
	for (i=1; i<=dim; i++)
		trannum[i] = gasdev(idum)/Y;

}

/*********************************************************************
 *   Returns the value ln(Gamma(xx)) for xx>0. Full accuracy is obtained
 *     for xx > 1. For 0<xx<1, the reflection formula can be used first:
 *
 *          Gamma(1-z) = pi/Gamma(z)/sin(pi*z) = pi*z/Gamma(1+z)/sin(pi*z)
 *          *********************************************************************/
double gammln(xx)
double xx;
{
        double x,tmp,ser;
        static double cof[6]={76.18009173,-86.50532033,24.01409822,
                -1.231739516,0.120858003e-2,-0.536382e-5};
        int j;

        x=xx-1.0;
        tmp=x+5.5;
        tmp -= (x+0.5)*log(tmp);
        ser=1.0;
        for (j=0;j<=5;j++) {
                x += 1.0;
                ser += cof[j]/x;
        }
        return -tmp+log(2.50662827465*ser);
}

#define PI 3.141592654

double poidev(xm,idum)
double xm;
long *idum;
{
        static double sq,alxm,g,oldm=(-1.0);
        double em,t,y;
        double ran2(),gammln();

        if (xm < 12.0) {
                if (xm != oldm) {
                        oldm=xm;
                        g=exp(-xm);
                }
                em = -1;
                t=1.0;
                do {
                        em += 1.0;
                        t *= ran2(idum);
                } while (t > g);
        } else {
                if (xm != oldm) {
                        oldm=xm;
                        sq=sqrt(2.0*xm);
                        alxm=log(xm);
                        g=xm*alxm-gammln(xm+1.0);
                }
                do {
                        do {
                                y=tan(PI*ran2(idum));
                                em=sq*y+xm;
                        } while (em < 0.0);
                        em=floor(em);
                        t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
                } while (ran2(idum) > t);
        }
        return em;
}

#undef PI


#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(idum)
long *idum;
{
        int j;
        long k;
        static long idum2=123456789;
        static long iy=0;
        static long iv[NTAB];
        double temp;

        if (*idum <= 0) {
                if (-(*idum) < 1) *idum=1;
                else *idum = -(*idum);
                idum2=(*idum);
                for (j=NTAB+7;j>=0;j--) {
                        k=(*idum)/IQ1;
                        *idum=IA1*(*idum-k*IQ1)-k*IR1;
                        if (*idum < 0) *idum += IM1;
                        if (j < NTAB) iv[j] = *idum;
                }
                iy=iv[0];
        }
        k=(*idum)/IQ1;
        *idum=IA1*(*idum-k*IQ1)-k*IR1;
        if (*idum < 0) *idum += IM1;
        k=idum2/IQ2;
        idum2=IA2*(idum2-k*IQ2)-k*IR2;
        if (idum2 < 0) idum2 += IM2;
        j=iy/NDIV;
        iy=iv[j]-idum2;
        iv[j] = *idum;
        if (iy < 1) iy += IMM1;
        if ((temp=AM*iy) > RNMX) return RNMX;
        else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
/* (C) Copr. 1986-92 Numerical Recipes Software 5.){2p491&.#@#Q0JS[. */

void sort(n,ra)
int n;
double ra[];
{
        int l,j,ir,i;
        double rra;

        l=(n >> 1)+1;
        ir=n;
        for (;;) {
                if (l > 1)
                        rra=ra[--l];
                else {
                        rra=ra[ir];
                        ra[ir]=ra[1];
                        if (--ir == 1) {
                                ra[1]=rra;
                                return;
                        }
                }
                i=l;
                j=l << 1;
                while (j <= ir) {
                        if (j < ir && ra[j] < ra[j+1]) ++j;
                        if (rra < ra[j]) {
                                ra[i]=ra[j];
                                j += (i=j);
                        }
                        else j=ir+1;
                }
                ra[i]=rra;
        }
}

void sort2(n,ra,rb)
int n;
double ra[],rb[];
{
        int l,j,ir,i;
        double rrb,rra;

        l=(n >> 1)+1;
        ir=n;
        for (;;) {
                if (l > 1) {
                        rra=ra[--l];
                        rrb=rb[l];
                } else {
                        rra=ra[ir];
                        rrb=rb[ir];
                        ra[ir]=ra[1];
                        rb[ir]=rb[1];
                        if (--ir == 1) {
                                ra[1]=rra;
                                rb[1]=rrb;
                                return;
                        }
                }
                i=l;
                j=l << 1;
                while (j <= ir) {
                        if (j < ir && ra[j] < ra[j+1]) ++j;
                        if (rra < ra[j]) {
                                ra[i]=ra[j];
                                rb[i]=rb[j];
                                j += (i=j);
                        }
                        else j=ir+1;
                }
                ra[i]=rra;
                rb[i]=rrb;
        }
}

void dmatrix_inv(double **A, double **invA, int n)
{
        double d, *col,**tempA;
        int i,j,*indx;

        col=dvector(1,n);
        indx=ivector(1,n);
        tempA=dmatrix(1,n,1,n);

        for (i=1;i<=n ;i++ )
                for (j=1;j<=n;j++)
                        tempA[i][j]=A[i][j];

        ludcmp(tempA,n,indx,&d);
        for (j=1;j<=n;j++)
        {
                for (i=1;i<=n;i++) col[i]=0.0;
                col[j]=1.0;
                lubksb(tempA,n,indx,col);
                for (i=1;i<=n;i++) invA[i][j]=col[i];
        }

        free_dvector(col,1,n);
        free_ivector(indx,1,n);
        free_dmatrix(tempA,1,n,1,n);

        return;
}

void choldc(double **a, int n, double *p)
/*
        given a positive-definite symmetric dmatrix a[1..n][1..n]. this routine constructs Cholesky
        decomposition. A=L * L'. On input, only the upper triangle of a need be given; it is not
        modified. The Cholesky factor L is returned in the lower triangle of a, except for its
        diagonal elements which are return in p[1..n]
*/
{
        int i,j,k,l,m;
        double sum;
        char temp;

        for (i=1;i<=n;i++) {
                for (j=i;j<=n;j++) {
                        for (sum=a[i][j],k=i-1;k>=1;k--) sum -= a[i][k]*a[j][k];
                        if (i == j)
			{
                                if (sum <= 0.0)
                                {
                                        printf("choldc failed\n");
	                                flag = 1;
                                }

                        p[i]=sqrt(sum);
                        } else a[j][i]=sum/p[i];
                }
        }
}


// Quasi-Newton Algorithm
#define ITMAX 1000
#define EPS 1.0e-10

void dfpmin(p,n,ftol,iter,fret,func,dfunc)
double p[],ftol,*fret,(*func)();
void (*dfunc)();
int n,*iter;
{
        int j,i,its;
        double fp,fae,fad,fac;
        double *xi,*g,*dg,*hdg,*dvector();
        double **hessin,**dmatrix();
        void linmin(),nrerror(),free_dmatrix(),free_dvector();

        hessin=dmatrix(1,n,1,n);
        xi=dvector(1,n);
        g=dvector(1,n);
        dg=dvector(1,n);
        hdg=dvector(1,n);
        fp=(*func)(p);
        (*dfunc)(p,g);
        for (i=1;i<=n;i++) {
                for (j=1;j<=n;j++) hessin[i][j]=0.0;
                hessin[i][i]=1.0;
                xi[i] = -g[i];
        }
        for (its=1;its<=ITMAX;its++) {

                *iter=its;
                linmin(p,xi,n,fret,func);
                if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
                        free_dvector(hdg,1,n);
                        free_dvector(dg,1,n);
                        free_dvector(g,1,n);
                        free_dvector(xi,1,n);
                        free_dmatrix(hessin,1,n,1,n);
                        return;
                }
                fp=(*fret);
                for (i=1;i<=n;i++) dg[i]=g[i];
                *fret=(*func)(p);
                (*dfunc)(p,g);
                for (i=1;i<=n;i++) dg[i]=g[i]-dg[i];
                for (i=1;i<=n;i++) {
                        hdg[i]=0.0;
                        for (j=1;j<=n;j++) hdg[i] += hessin[i][j]*dg[j];
                }
                fac=fae=0.0;
                for (i=1;i<=n;i++) {
                        fac += dg[i]*xi[i];
                        fae += dg[i]*hdg[i];
                }
                fac=1.0/fac;
                fad=1.0/fae;
                for (i=1;i<=n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
                for (i=1;i<=n;i++)
                        for (j=1;j<=n;j++)
                                hessin[i][j] += fac*xi[i]*xi[j]
                                        -fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
                for (i=1;i<=n;i++) {
                        xi[i]=0.0;
                        for (j=1;j<=n;j++) xi[i] -= hessin[i][j]*g[j];
                }
        }
        /*
        nrerror("Too many iterations in DFPMIN");
        */
}

#undef ITMAX
#undef EPS
#define TOL 2.0e-4

int ncom=0;     /* defining declarations */
double *pcom=0,*xicom=0,(*nrfunc)();

void linmin(p,xi,n,fret,func)
double p[],xi[],*fret,(*func)();
int n;
{
        int i,j;
        double xx,xmin,fx,fb,fa,bx,ax;
        double brent(),f1dim(),*dvector();
        void mnbrak(),free_dvector();
        double maxxi;

        ncom=n;
        pcom=dvector(1,n);
        xicom=dvector(1,n);
        nrfunc=func;
        for (j=1;j<=n;j++) {
                pcom[j]=p[j];
                xicom[j]=xi[j];
        }
        ax=0.0;
        xx=1.0e-3;
        bx=2.0;
        mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
        *fret=brent(ax,xx,bx,f1dim,TOL,&xmin);

        maxxi = fabs(xi[1]);
        for (j=2; j<=n; j++)
                if (fabs(xi[j])>maxxi)
                        maxxi = fabs(xi[j]);

        for (j=1;j<=n;j++) {
                xi[j] *= xmin;
                p[j] += xi[j];
        }

        free_dvector(xicom,1,n);
        free_dvector(pcom,1,n);
}

#undef TOL

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define MAX(a,b) ((a) > (b) ? (a) : (b))
//#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

void mnbrak(ax,bx,cx,fa,fb,fc,func)
double *ax,*bx,*cx,*fa,*fb,*fc;
double (*func)();       /* ANSI: double (*func)(double); */
{
        double ulim,u,r,q,fu,dum;

        *fa=(*func)(*ax);
        *fb=(*func)(*bx);
        if (*fb > *fa) {
                SHFT(dum,*ax,*bx,dum)
                SHFT(dum,*fb,*fa,dum)
        }
        *cx=(*bx)+GOLD*(*bx-*ax);
        *fc=(*func)(*cx);
        while (*fb > *fc) {
                r=(*bx-*ax)*(*fb-*fc);
                q=(*bx-*cx)*(*fb-*fa);
                u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
                        (2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
                ulim=(*bx)+GLIMIT*(*cx-*bx);
                if ((*bx-u)*(u-*cx) > 0.0) {
                        fu=(*func)(u);
                        if (fu < *fc) {
                                *ax=(*bx);
                                *bx=u;
                                *fa=(*fb);
                                *fb=fu;
                                return;
                        } else if (fu > *fb) {
                                *cx=u;
                                *fc=fu;
                                return;
                        }
                        u=(*cx)+GOLD*(*cx-*bx);
                        fu=(*func)(u);
                } else if ((*cx-u)*(u-ulim) > 0.0) {
                        fu=(*func)(u);
                        if (fu < *fc) {
                                SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
                                SHFT(*fb,*fc,fu,(*func)(u))
                        }
                } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
                        u=ulim;
                        fu=(*func)(u);
                } else {
                        u=(*cx)+GOLD*(*cx-*bx);
                        fu=(*func)(u);
                }
                SHFT(*ax,*bx,*cx,u)
                SHFT(*fa,*fb,*fc,fu)
        }
}

#undef GOLD
#undef GLIMIT
#undef TINY
#undef MAX
//#undef SIGN
#undef SHFT
extern int ncom;        /* defined in LINMIN */
extern double *pcom,*xicom,(*nrfunc)();

double f1dim(x)
double x;
{
        int j;
        double f,*xt,*dvector();
        void free_dvector();

        xt=dvector(1,ncom);
        for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];

        f=(*nrfunc)(xt);
        free_dvector(xt,1,ncom);
        return f;
}

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
//#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))^M
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

double brent(ax,bx,cx,f,tol,xmin)
double ax,bx,cx,tol,*xmin;
double (*f)();  /* ANSI: double (*f)(double); */
{
        int iter;
        double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
        double e=0.0;
        void nrerror();

        a=((ax < cx) ? ax : cx);
        b=((ax > cx) ? ax : cx);
        x=w=v=bx;
        fw=fv=fx=(*f)(x); //printf("fx=%f\n", fx);
        for (iter=1;iter<=ITMAX;iter++) {
                xm=0.5*(a+b);
                tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
                if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
                        *xmin=x;
                        return fx;
                }
                if (fabs(e) > tol1) {
                        r=(x-w)*(fx-fv);
                        q=(x-v)*(fx-fw);
                        p=(x-v)*q-(x-w)*r;
                        q=2.0*(q-r);
                        if (q > 0.0) p = -p;
                        q=fabs(q);
                        etemp=e;
                        e=d;
                        if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
                                d=CGOLD*(e=(x >= xm ? a-x : b-x));
                        else {
                                d=p/q;
                                u=x+d;
                                if (u-a < tol2 || b-u < tol2)
                                        d=SIGN(tol1,xm-x);
                        }
                } else {
                        d=CGOLD*(e=(x >= xm ? a-x : b-x));
                }
                u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
                fu=(*f)(u);
                if (fu <= fx) {
                        if (u >= x) a=x; else b=x;
                        SHFT(v,w,x,u)
                        SHFT(fv,fw,fx,fu)
                } else {
                        if (u < x) a=u; else b=u;
                        if (fu <= fw || w == x) {
                                v=w;
                                w=u;
                                fv=fw;
                                fw=fu;
                        } else if (fu <= fv || v == x || v == w) {
                                v=u;
                                fv=fu;
                        }
                }
        }
        nrerror("Too many iterations in BRENT");
        *xmin=x;
        return fx;
}

#undef ITMAX
#undef CGOLD
#undef ZEPS
//#undef SIGN
