//#include <mnewt.h>
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <Rmath.h>
#include <stdio.h>


/* pre-difine something we need */
#define NR_END 1
#define FREE_ARG char*

#define TINY 1.0e-20
void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

/* some data tpye */
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

   
double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
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

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
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
 

// pre-define functions of mybspline_c.c
void ludcmp(double **a, int n, int *indx, double *d)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;

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


double distance(float *x1, float *x2, int p)
{
    int m;
    double res=0;
    for(m=0;m<p;m++)
        res=res+pow(x1[m+1]-x2[m+1],2);
    res=sqrt(res);
    return(res);
}




// main body of mybspline_c.c
SEXP compbmymmc(SEXP X, SEXP Y, SEXP TAU, SEXP THETA, SEXP N, SEXP PM, SEXP KN, SEXP MAXITER, SEXP TOL, SEXP EPSILON)
{
    // X: enlarged predictor matrix
    // Y: enlarged response vector
    // TAU: enlarged tau vector and these quantile points are the ones that we are intrersted
    // THETA: starting point of the coef
    // N: number of original observations
    // PM: dimensions of enlarged predictors, equals to p*m where m is the number of basis functions
    // KN: number of TAUs(in our paper ,KN+1 is the number of TAUs)
    // MAXITER: max number of iterations
    // TOL: tolerance of convergence
    // EPSILON: parameter used in MM
    
    int i,j,m;   // index for loops
    int l=0;   // number of MM iterations
    int n=INTEGER(N)[0];
    int pm=INTEGER(PM)[0];
    int kn=INTEGER(KN)[0];
    int maxiter=INTEGER(MAXITER)[0];
    double tol=REAL(TOL)[0];
    double epsilon=REAL(EPSILON)[0];
    double **x, *y, *tau, *theta;
    double err=100.0;   // store the sum of absolute difference between theta and theta.old
    double *r;   // store the residual vector
    double *theta_old;    // store the theta of last loop, so we can compare the difference
    double *denominator;   // store the vector 1/(epsilon + abs(r_ik) )
    double *v;    // store the vector 4*\tau - 2
    // In this algorithm, we update theta by solving the linear equations:
    // (Z^T A Z) theta = Z^T(A Y + 0.5 B)
    // where Z is matrix x and Y is vector y in this program
    // A is a diagnoal matrix, with its diagnoal entries stored in denominator in this program
    // B is a vector, which is v in this program
    double **zaz, *za;    // zaz corresbons to (Z^T A Z) and za corresponds to Z^T(A Y + 0.5 B)
    double d;   // parameter needed by ludcomp
    int *index;   // parameter neeeded by ludcomp
    
    SEXP RES=PROTECT(allocVector(REALSXP,pm));   // communicate resutls with R
    
    // initialize the forms
    x=dmatrix(1,n*kn,1,pm);
    y=dvector(1,n*kn);
    tau=dvector(1,n*kn);
    theta=dvector(1,pm);
    theta_old=dvector(1,pm);
    r=dvector(1,n*kn);
    v=dvector(1,n*kn);
    denominator=dvector(1,n*kn);
    index=ivector(1,pm);
    zaz=dmatrix(1,pm,1,pm);
    za=dvector(1,pm);
    // initialize the values
    for(i=0;i<n*kn;i++)
        for(j=0;j<pm;j++)
            x[i+1][j+1]=REAL(X)[j*n*kn+i];
    for(i=0;i<n*kn;i++)   
    {
        y[i+1]=REAL(Y)[i];
        tau[i+1]=REAL(TAU)[i];
        v[i+1]=4.0*tau[i+1]-2.0;    // this value is unchanged during the algorithm
    }   
    for(i=0;i<pm;i++)
    {
        theta[i+1]=REAL(THETA)[i];
        theta_old[i+1]=theta[i+1];
    }        

    
    // MM algorithm
    while(err>tol && l<maxiter)
    {
        // save the last theta 
        for(i=0;i<pm;i++)   theta_old[i+1]=theta[i+1];
        
        // compute r, denominator and v
        for(i=0;i<n*kn;i++)
        {
            // compute residual vector
            r[i+1]=y[i+1];
            for(j=0;j<pm;j++)   r[i+1]=r[i+1]-x[i+1][j+1]*theta[j+1];
            
            // compute denominator vector
            denominator[i+1]=1.0/( epsilon + fabs(r[i+1]) );                    
        }
        
        // compute za ( Z^T(A Y + 0.5 B) in the draft )
        for(i=0;i<pm;i++)
        {
            za[i+1]=0;
            for(j=0;j<n*kn;j++)   za[i+1]=za[i+1] + x[j+1][i+1]*(denominator[j+1]*y[j+1]+0.5*v[j+1]);
        }
        
        // compute zaz ( (Z^T A Z) )
        for(i=0;i<pm;i++)
            for(j=0;j<pm;j++)
            {
                zaz[i+1][j+1]=0;
                for(m=0;m<n*kn;m++)   zaz[i+1][j+1]=zaz[i+1][j+1] + x[m+1][i+1]*denominator[m+1]*x[m+1][j+1];
            }
        
        // compute delta
        ludcmp(zaz,pm,index,&d);   
        lubksb(zaz,pm,index,za);   // after this , the value of za has become the counterpart of theta
        
        // update theta, err and l
        err=0;
        for(i=0;i<pm;i++)
        {
            theta[i+1]=za[i+1];
            err=err+fabs(theta[i+1]-theta_old[i+1]);
        }            
        l=l+1;     
    }
     printf("%d\t %.10f\t %.10f\t %.10f\n",l,err,tol,epsilon);
    
    // prepare the result
    for(i=0;i<pm;i++)   REAL(RES)[i]=theta[i+1];
    // free the allocation
    free_dmatrix(x,1,n*kn,1,pm);
    free_dvector(y,1,n*kn);
    free_dvector(tau,1,n*kn);
    free_dvector(theta,1,pm);
    free_dvector(theta_old,1,pm);
    free_dvector(r,1,n*kn);
    free_dmatrix(zaz,1,pm,1,pm);
    free_dvector(za,1,pm);
    free_dvector(v,1,n*kn);
    free_dvector(denominator,1,n*kn);
    free_ivector(index,1,pm);
    UNPROTECT(1);
    // return the result
    return(RES);
}


SEXP compbmmc(SEXP X, SEXP Y, SEXP TAU, SEXP THETA, SEXP N, SEXP PM, SEXP KN, SEXP MAXITER, SEXP TOL, SEXP EPSILON)
{
    // X: enlarged predictor matrix
    // Y: enlarged response vector
    // TAU: enlarged tau vector and these quantile points are the ones that we are intrersted
    // THETA: starting point of the coef
    // N: number of original observations
    // PM: dimensions of enlarged predictors, equals to p*m where m is the number of basis functions
    // KN: number of TAUs(in our paper ,KN+1 is the number of TAUs)
    // MAXITER: max number of iterations
    // TOL: tolerance of convergence
    // EPSILON: parameter used in MM
    
    int i,j,m;   // index for loops
    int l=0;   // number of MM iterations
    int n=INTEGER(N)[0];
    int pm=INTEGER(PM)[0];
    int kn=INTEGER(KN)[0];
    int maxiter=INTEGER(MAXITER)[0];
    double tol=REAL(TOL)[0];
    double epsilon=REAL(EPSILON)[0];
    double **x, *y, *tau, *theta;
    double err=100.0;   // store the sum of absolute difference between theta and theta.old
    double *r;   // store the residual vector
    //float *delta;   // the difference between theta and theta_old 
    // delta is the solution for  X^TWX \cdot delta = X^T\cdot V
    // by design of the lubksb, the value of delta will be stored in xv
    // in this programe, theta.new = theta.old - delta = theta.old - xv
    double **xwx;   // left matrix of the delta equation, it is a pm * pm matrix
    double *xv;   // right vector of the delta equation, it is a pm*1 vector
    // float *theta_old;   // store the old value of theta
    double *v;   // store the vector V
    double *denominator;   // store the vector 1/(epsilon + abs(r_ik) )
    double d;   // parameter needed by ludcomp
    int *index;   // parameter neeeded by ludcomp
    
    SEXP RES=PROTECT(allocVector(REALSXP,pm));   // communicate resutls with R
    
    // initialize the forms
    x=dmatrix(1,n*kn,1,pm);
    y=dvector(1,n*kn);
    tau=dvector(1,n*kn);
    theta=dvector(1,pm);
    
    r=dvector(1,n*kn);
    // delta=vector(1,pm);
    xwx=dmatrix(1,pm,1,pm);
    xv=dvector(1,pm);
    // theta_old=vector(1,pm);
    v=dvector(1,n*kn);
    denominator=dvector(1,n*kn);
    index=ivector(1,pm);
    
    // initialize the values
    for(i=0;i<n*kn;i++)
        for(j=0;j<pm;j++)
            x[i+1][j+1]=REAL(X)[j*n*kn+i];
    for(i=0;i<n*kn;i++)   
    {
        y[i+1]=REAL(Y)[i];
        tau[i+1]=REAL(TAU)[i];
    }   
    for(i=0;i<pm;i++)   theta[i+1]=REAL(THETA)[i];

    
    // MM algorithm
    while(err>tol && l<maxiter)
    {
        // save the last theta 
        // for(i=0;i<pm;i++)   theta_old[i+1]=theta[i+1];
        
        // compute r, denominator and v
        for(i=0;i<n*kn;i++)
        {
            // compute residual vector
            r[i+1]=y[i+1];
            for(j=0;j<pm;j++)   r[i+1]=r[i+1]-x[i+1][j+1]*theta[j+1];
            
            // compute denominator vector
            denominator[i+1]=1.0/( epsilon + fabs(r[i+1]) );
            
            // compute V vector
            v[i+1]= 1.0-2.0*tau[i+1]-r[i+1]*denominator[i+1];         
        }
        
        // compute xv
        for(i=0;i<pm;i++)
        {
            xv[i+1]=0;
            for(j=0;j<n*kn;j++)   xv[i+1]=xv[i+1]+x[j+1][i+1]*v[j+1];
        }
        
        // compute xwx
        for(i=0;i<pm;i++)
            for(j=0;j<pm;j++)
            {
                xwx[i+1][j+1]=0;
                for(m=0;m<n*kn;m++)   xwx[i+1][j+1]=xwx[i+1][j+1] + x[m+1][i+1]*denominator[m+1]*x[m+1][j+1];
            }
        
        // compute delta
        ludcmp(xwx,pm,index,&d);   
        lubksb(xwx,pm,index,xv);   // after this , the value of xv has become the counterpart of delta
        
        // update theta, err and l
        err=0;
        for(i=0;i<pm;i++)
        {
            theta[i+1]=theta[i+1]-xv[i+1];
            err=err+fabs(xv[i+1]);
        }            
        l=l+1;     
    }
     printf("%d\t %.10f\t %.10f\t %.10f\n",l,err,tol,epsilon);
    
    // prepare the result
    for(i=0;i<pm;i++)   REAL(RES)[i]=theta[i+1];
    // free the allocation
    free_dmatrix(x,1,n*kn,1,pm);
    free_dvector(y,1,n*kn);
    free_dvector(tau,1,n*kn);
    free_dvector(theta,1,pm);
    free_dvector(r,1,n*kn);
    free_dmatrix(xwx,1,pm,1,pm);
    free_dvector(xv,1,pm);
    free_dvector(v,1,n*kn);
    free_dvector(denominator,1,n*kn);
    free_ivector(index,1,pm);
    UNPROTECT(1);
    // return the result
    return(RES);
}


SEXP unc(SEXP X, SEXP PHI, SEXP N, SEXP P)  // compute LOF test statistic
{
    // X: orginal observations( only predictors )
    // PHI: results from function Phi
    // N: number of observations, which is nrow(x) and nrow(phi)
    // P: dimensions of predictors, which is ncol(x)
    
    int i,j,k;
    int n=INTEGER(N)[0];
    int p=INTEGER(P)[0];
    float **x, *phi;
    float *x1, *x2;  // store 2 predictors in order to compute the distance between them
    double dis;   // store the distance of 2 observations
    double res;   // store the test statistic
    
    SEXP RES=PROTECT(allocVector(REALSXP,1));   // comunicate results with R
    
    // initialize the forms
    x=matrix(1,n,1,p);
    phi=vector(1,n);
    x1=vector(1,p);
    x2=vector(1,p);
    
    // initialize the values 
    for(i=0;i<n;i++) phi[i+1]=REAL(PHI)[i];
    for(i=0;i<n;i++)
        for(j=0;j<p;j++)
            x[i+1][j+1]=REAL(X)[j*n+i];
    
    // compute the test statistic
    res=0;
    for(i=0;i<n;i++)
        for(j=i+1;j<n;j++)
        {
            for(k=0;k<p;k++)
            {
                x1[k+1]=x[i+1][k+1];
                x2[k+1]=x[j+1][k+1];
            }
            dis=distance(x1,x2,p);
            res=res+dis*phi[i+1]*phi[j+1];
        }
    res=res*2/n/(n-1);
    res=fabs(res);   // compute the absolute value
    
    // prepare the result
    REAL(RES)[0]=res;
    // free the allocations
    free_matrix(x,1,n,1,p);
    free_vector(phi,1,n);
    free_vector(x1,1,p);
    free_vector(x2,1,p);
    UNPROTECT(1);
    // return the result
    return(RES);
        
}



SEXP unbc(SEXP X, SEXP X_B, SEXP PHI, SEXP PHI_B, SEXP N, SEXP P)   // compute bootstrap LOF test statistic
{
    // X: original observation( only predictors )
    // X.B: bootstrap observation
    // PHI: results of function Phi with respect to X
    // PHI.B: results of function Phi with respect to X.B
    // N: number of observations
    // P: dimensions of predictor
    
    int i,j,k,l,m;   // index used in for loops
    int n=INTEGER(N)[0];
    int p=INTEGER(P)[0];
    float **x, **x_b, *phi, *phi_b;
    float *x1,*x2;   // store 2 predictors in order to compute the distance between them
    float *ksum;   // store the \sum\limits_{k=1}^n between x and x.b
    double dis=0;   // store the disntance between 2 predictors
    double res=0, res1=0, res2=0;   // store the test statistics 
    double tmp1=0, tmp2=0, tmp3=0;   // store some temperature values during the computation    
    SEXP RES=PROTECT(allocVector(REALSXP,1));   // communicate results with R
    
    // initialize some forms
    x=matrix(1,n,1,p);
    x_b=matrix(1,n,1,p);
    phi=vector(1,n);
    phi_b=vector(1,n);
    ksum=vector(1,n);
    x1=vector(1,p);
    x2=vector(1,p);
    // initialize some values
    for(i=0;i<n;i++)
    {
        phi[i+1]=REAL(PHI)[i];
        phi_b[i+1]=REAL(PHI_B)[i];
    }  
    for(i=0;i<n;i++)        
        for(j=0;j<p;j++)
        {
            x[i+1][j+1]=REAL(X)[j*n+i];
            x_b[i+1][j+1]=REAL(X_B)[j*n+i];
        }
    
    // compute test statistic
    // compute the first part of the statistic
    res1=0;  
    for(i=0;i<n;i++)
        for(j=i+1;j<n;j++)
        {
            for(m=0;m<p;m++)
            {
                x1[m+1]=x_b[i+1][m+1];
                x2[m+1]=x_b[j+1][m+1];
            }
            dis=distance(x1,x2,p);
            
            // update part 1 in every loop
            res1=res1+phi_b[i+1]*phi_b[j+1]*dis;
        }
    res1=res1*2/n/(n-1);
    
    // compute the second part of the statistic
    res2=0;
    // since the third part (\sum\limits_{k,l=1}^n) is the same in every for loop of i and j
    // compute it here first
    tmp3=0;   // compute the third \sum\limits_{k,l=1}^n
    for(k=0;k<n;k++)
        for(l=0;l<n;l++)
        {
            for(m=0;m<p;m++)
            {
                x1[m+1]=x[k+1][m+1];
                x2[m+1]=x[l+1][m+1];
            }
            dis=distance(x1,x2,p);  
            tmp3=tmp3+phi[k+1]*phi[l+1]*dis;                 
        } 
    tmp3=tmp3/n;
    
    // since \sum\limits_{k=1}^n between x and x.b will be used all the time 
    // compute it here and store in ksum
    for(i=0;i<n;i++)
    {
        ksum[i+1]=0;
        for(k=0;k<n;k++)
        {
            for(m=0;m<p;m++)
            {
                x1[m+1]=x_b[i+1][m+1];
                x2[m+1]=x[k+1][m+1];
            }
            dis=distance(x1,x2,p);
            ksum[i+1]=ksum[i+1]+phi[k+1]*dis;
        }
    }
    
    for(i=0;i<n;i++)
        for(j=i+1;j<n;j++)
        {
            // compute the first \sum\limits_{k=1}^n with x_i^\star and x_k
            tmp1=phi_b[i+1]*ksum[i+1];
            
            // compute the second \sum\limist_{k=1}^n with x_k and x_j^\star
            tmp2=phi_b[j+1]*ksum[j+1];
            
            // update part 2 in every loop
            res2=res2+tmp1+tmp2-tmp3;            
        }
    res2=res2*2/n/n/(n-1);
    
    res=res1-res2;
    res=fabs(res);   //compute the absolute value
    
    // prepare the result
    REAL(RES)[0]=res;
    // free the allocations
    free_matrix(x,1,n,1,p);
    free_matrix(x_b,1,n,1,p);
    free_vector(phi,1,n);
    free_vector(phi_b,1,n);
    free_vector(ksum,1,n);
    free_vector(x1,1,p);
    free_vector(x2,1,p);
    UNPROTECT(1);
    // return the result
    return(RES);
}

#undef TINY