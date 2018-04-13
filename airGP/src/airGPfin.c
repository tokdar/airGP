#include "R.h"
#include "Rmath.h"
#include "R_ext/Applic.h"
#include "time.h"

// a library of utility functions

double *vect(int n);
void vect_set(double *x, int n, double val);
void vect_import(double *x, int n, double *y);
int *ivect(int n);
double **mymatrix(int nr, int nc);
void Rprintvec(char *a, char *format, double *x, int n);
void Rprintmat(char *a, char *format, double **x, int m, int n, int flip);
void Rprintveci(char *a, char *format, int *x, int n);
double sumsquares(double *x, int n);
double inprod(double *x, double *y, int n);
double inprod_subset(double *x, double *y, int n, int *ix);
double rnormtrunc(double mu, double sigma, double lo, double hi);
double rgammatrunc(double shape, double scale, double lo, double hi);
double rbetatrunc(double shape1, double shape2, double lo, double hi);
double vmax(double *x, int n);
double vmin(double *x, int n);
double logsum(double *lx, int n);
double sum(double *x, int n);
double sum_subset(double *x, int n, int *ix);
double sumf(double *x, int n, double f(double));
double sumf_subset(double *x, int n, int *ix, double f(double));
int isum(int *x, int n);
double average(double *x, int n);
double variance(double *x, int n);
double sigmoid_lin(double x);
double adjustmentFn(double x);
int rdraw(int n, double *lprob, int inlog);
double rdraw_trapez(double *f, double *t, int n, double *Fn);
double rposnorm(double mu, double sigma);
double dcauchy(double x, double mu, double sigma, int islog);
void locator_string(int *ix, int n, char *a);
void string_locator(char *a, int *ix, int n);
int any(int *x, int n);
void q_sort(double *xinput, int *xpos, int xlen, int start);
void quick_sort(double *xinput, int *xpos, int xlen, int start);
double matern(double x, double phi, int kappa);
void set_lower_tri_zero(double **A, int n, int m);
void spchol_st(double **R, int N, double tol, int *pivot, int *rank, int max_rank, double gpcov(int, int, double*, int*), double *covpar, int *include, int max_scan, double *d, int dopivoting, int dostopping, int padzero, double *lpen);
void spchol_agp(double **R, int N, double tol, int *pivot, int *rank, int max_rank, double gpcov(int, int, double*, int*), int ncomp, int npar, int nvar, double *covpar, int *include, int max_scan, double *d, int dopivoting, int dostopping, int padzero, double *lpen, int *active);
void chol(double **R, int N, double tol, int *pivot, int *rank, int max_rank, double *d, double **A, int dopivoting, int padzero);
void trisolve(double **R, int m, double *b, double *x, int transpose);
void triprod(double **R, int m, int n, double *x, double *b, int transpose);
double stabler(double x, double threshold);

double * vect(int n){
    return (double *)R_alloc(n, sizeof(double));
}
void vect_set(double *x, int n, double val){
    int i;
    for(i = 0; i < n; i++) x[i] = val;
}
void vect_import(double *x, int n, double *y){
    int i;
    for(i = 0; i < n; i++) x[i] = y[i];
}
int * ivect(int n){
    return (int *)R_alloc(n, sizeof(int));
}
double ** mymatrix(int nr, int nc){
    int   i;
    double **m;
    m = (double **) R_alloc(nr, sizeof(double *));
    for (i = 0; i < nr; i++)
        m[i] = (double *) R_alloc(nc, sizeof(double));
    return m;
}
void Rprintvec(char *a, char *format, double *x, int n){
    int i;
    Rprintf("%s", a);
    for(i = 0; i < n; i++)
        Rprintf(format, x[i]);
    Rprintf("\n");
}
void Rprintmat(char *a, char *format, double **x, int m, int n, int flip){
    int i, j;
    Rprintf("%s\n", a);
    for(i = 0; i < m; i++){
        for(j = 0; j < n; j++)
            Rprintf(format, x[i][j]);
        Rprintf("\n");
    }
}
void Rprintveci(char *a, char *format, int *x, int n){
    int i;
    Rprintf("%s", a);
    for(i = 0; i < n; i++)
        Rprintf(format, x[i]);
    Rprintf("\n");
}
double inprod(double *x, double *y, int n){
    double ip = 0.0;
    int i;
    for(i = 0; i < n; i++)
        ip += x[i] * y[i];
    return ip;
}
double inprod_subset(double *x, double *y, int n, int *ix){
    double ip = 0.0;
    int i;
    for(i = 0; i < n; i++)
        if(ix[i])
            ip += x[i] * y[i];
    return ip;
}
double sumsquares(double *x, int n){
    return inprod(x, x, n);
}
double rnormtrunc(double mu, double sigma, double lo, double hi){
    double u = runif(0.0, 1.0);
    double p = u * pnorm(hi, mu, sigma, 1, 0) + (1.0 - u) * pnorm(lo, mu, sigma, 1, 0);
    if(p <= 0.0) p = 1.0e-10;
    if(p >= 1.0) p = 1.0 - 1.0e-10;
    return qnorm(p, mu, sigma, 1, 0);
}
double rgammatrunc(double shape, double scale, double lo, double hi){
    double u = runif(0.0, 1.0);
    double p = u * pgamma(hi, shape, scale, 1, 0) + (1.0 - u) * pgamma(lo, shape, scale, 1, 0);
    if(p <= 0.0) p = 1.0e-10;
    if(p >= 1.0) p = 1.0 - 1.0e-10;
    return qgamma(p, shape, scale, 1, 0);
}
double rbetatrunc(double shape1, double shape2, double lo, double hi){
    double u = runif(0.0, 1.0);
    double p = u * pbeta(hi, shape1, shape2, 1, 0) + (1.0 - u) * pbeta(lo, shape1, shape2, 1, 0);
    if(p <= 0.0) p = 1.0e-10;
    if(p >= 1.0) p = 1.0 - 1.0e-10;
    return qbeta(p, shape1, shape2, 1, 0);
}
double vmax(double *x, int n){
    int i;
    double xmax = x[0];
    for(i = 1; i < n; i++) if(x[i] > xmax) xmax = x[i];
    return xmax;
}
double vmin(double *x, int n){
    int i;
    double xmin = x[0];
    for(i = 1; i < n; i++) if(x[i] < xmin) xmin = x[i];
    return xmin;
}
double logsum(double *lx, int n){
    double lxmax = vmax(lx, n), a = 0.0;
    int i;
    for(i = 0; i < n; i++) a += exp(lx[i] - lxmax);
    return lxmax + log(a);
}
double sum(double *x, int n){
    double a = 0.0;
    int i;
    for(i = 0; i < n; i++) a += x[i];
    return a;
}
double sum_subset(double *x, int n, int *ix){
    int i;
    double a = 0.0;
    for(i = 0; i < n; i++) if(ix[i]) a += x[i];
    return a;
}
double sumf(double *x, int n, double f(double)){
    int i;
    double a = 0.0;
    for(i = 0; i < n; i++) a += f(x[i]);
    return a;
}
double sumf_subset(double *x, int n, int *ix, double f(double)){
    int i;
    double a = 0.0;
    for(i = 0; i < n; i++) if(ix[i]) a += f(x[i]);
    return a;
}
int isum(int *x, int n){
    int i, val = 0;
    for(i = 0; i < n; i++) val += x[i];
    return val;
}
double average(double *x, int n){
    return sum(x, n) / ((double)n);
}
double variance(double *x, int n){
    double xmean = average(x, n);
    return sumsquares(x, n) / (double)n - xmean*xmean;
}
double sigmoid_lin(double x){
    return x / (1.0 + x);
}
double adjustmentFn(double x){
    return x;
}
int rdraw(int n, double *prob, int inlog){
    double psum, u = runif(0.0, 1.0), cprob;
    int j = 0;
    
    if(inlog){
        psum = logsum(prob, n);
        cprob = exp(prob[0] - psum);
        while(u > cprob && j < n - 1){
            j++;
            cprob += exp(prob[j] - psum);
        }
    } else {
        psum = sum(prob, n);
        cprob = prob[0] / psum;
        while(u > cprob && j < n - 1){
            j++;
            if(prob[j] > 0.0) cprob += prob[j] / psum;
        }
    }
    return j;
}
double rdraw_trapez(double *f, double *t, int n, double *Fn){
    int i;
    double v = runif(0.0, 1.0), Ft = 0.0;
    for(i = 1; i < n; i++) Ft += 0.5*(f[i-1]+f[i])*(t[i]-t[i-1]);
    Fn[0] = Ft;
    v *= Ft;
    i = 1;
    double incr = 0.5*(f[i-1]+f[i])*(t[i]-t[i-1]);
    while(v > incr){
        v -= incr;
        i++;
        incr = 0.5*(f[i-1]+f[i])*(t[i]-t[i-1]);
    }
    double r = 2.0*v/(f[i-1] + sqrt(f[i-1]*f[i-1] + 2.0*(f[i] - f[i-1])*v));
    return t[i-1] + (t[i]-t[i-1])*r;
}
double rposnorm(double mu, double sigma){
    double u = runif(0.0, 1.0);
    return qnorm(u + (1 - u) * pnorm(0.0, mu, sigma, 1, 0), mu, sigma, 1, 0);
}
double dcauchy(double x, double mu, double sigma, int islog){
    double z = (x - mu)/sigma;
    double lf = -log(PI) - log(sigma) - log1p(z*z);
    if(!islog) lf = exp(lf);
    return lf;
}
void locator_string(int *ix, int n, char *a){
    const char *fmt[2];	fmt[0] = "%d";	fmt[1] = ".%d";
    int i, skip = 0;
    for(i = 0; i < n; i++){
        if(ix[i]){
            sprintf(a + skip, fmt[skip > 0], i + 1);
            skip = strlen(a);
        }
    }
}
void string_locator(char *a, int *ix, int n){
    char seps[] = ".";
    char* token;
    int var;
    
    token = strtok (a, seps);
    while (token != NULL)
    {
        sscanf (token, "%d", &var);
        ix[var - 1] = 1;
        token = strtok (NULL, seps);
    }
    
}
int any(int *x, int n){
    int i = 0;
    while(!x[i] && i < n) i++;
    return (i < n);
}
double stabler(double x, double threshold){
    if(x < threshold) x = threshold;
    return x;
}
void q_sort(double *xinput, int *xpos, int xlen, int start){
    int j;
    for(j = 0; j < xlen; j++) xpos[j] = j;
    quick_sort(xinput, xpos, xlen, start);
}

void quick_sort(double *xinput, int *xpos, int xlen, int start){
    
    int i, j, posdummy;
    
    if(xlen > 1){                                    // Do nothing if x is of length 1 or less
        if(start){
            i = floor(xlen * runif(0.0,1.0));        // Randomly select a splitter
            posdummy = xpos[0];                      // Swap to bring the splitter
            xpos[0] = xpos[i];                       // to the first position
            xpos[i] = posdummy;
        }
        
        i = 0;                                      // i stores splitter position
        j = 0;                                      // j stores scan range
        
        while(++j < xlen){                          // Increase scan length
            if(xinput[xpos[j]] < xinput[xpos[i]]){  // If the next scanned is smaller than the splitter then
                posdummy = xpos[i];                 // swap splitter to right, splitter's next right neighbor
                xpos[i] = xpos[j];                  // to the scanned and the scanned to the splitter
                i++;
                xpos[j] = xpos[i];
                xpos[i] = posdummy;
            }
        }                                           // This will put the splitter in its right place
        
        
        quick_sort(xinput, xpos, i, start);                 // Perform quick sort separately on the two subsequences //
        quick_sort(xinput, xpos + i + 1, xlen - i - 1, start); // formed by the splitter                                //
    }
}

double matern(double x, double phi, int kappa){
    /*Returns the Matern function for x, phi and kappa.*/
    
    /* Variables */
    double ans, cte;
    double uphi=x/phi;
    
    /* Matern */
    
    if (uphi==0) return 1;
    else{
        if (kappa==0.5)
            ans = exp(-uphi);
        else {
            cte = R_pow(2, (-(kappa-1)))/gammafn(kappa);
            ans = cte * R_pow(uphi, kappa) * bessel_k(uphi,kappa,1);
        }
    }
    
    return ans;
}

// Sparse Cholesky factorization with Pivoting
void set_lower_tri_zero(double **A, int n, int m ){
    int i, j;
    for(i = 0; i < n; i++)
        for(j = i + 1; j < m; j++)
            A[j][i] = 0.0;
}

void spchol_st(double **R, int N, double tol, int *pivot, int *rank, int max_rank,
               double gpcov(int, int, double*, int*), double *covpar, int *include,
               int max_scan, double *d, int dopivoting, int dostopping, int padzero, double *lpen){
    
    
    // Sparse cholesky factorization with pivoting and diagonal augmentation
    // for stationary covariance functions.
    // Accepts an empty matrix R and a function gpcov(i, j, ...) that is called
    // to compute the (i,j)-th element of the original covariance matrix
    
    
    
    set_lower_tri_zero(R, N, max_rank);
    
    int i, a, l;
    double u, b;
    int randFirst = 1;
    for(i = 0; i < N; i++) pivot[i] = i;
    if(dopivoting) q_sort(lpen, pivot, max_scan, randFirst);
    
    for(i = 0; i < N; i++) d[i] = gpcov(pivot[i], pivot[i], covpar, include);
    
    int k = 0, max_diag;
    for(max_diag = k, i = k + 1; i < max_scan; i++)
        if(d[i] > d[max_diag])
            max_diag = i;
    tol *= d[max_diag];
    int flag = (k < max_rank);
    if(dostopping)
        flag = (d[max_diag] > tol);
    
    while(flag){
        if(dopivoting){
            if(max_diag > k){
                a = pivot[k];
                pivot[k] = pivot[max_diag];
                pivot[max_diag] = a;
                
                b = d[k];
                d[k] = d[max_diag];
                d[max_diag] = b;
                
                for(i = 0; i < k; i++){
                    b = R[i][k];
                    R[i][k] = R[i][max_diag];
                    R[i][max_diag] = b;
                }
            }
        }
        
        R[k][k] = sqrt(d[k]);
        
        for(i = k + 1; i < N; i++){
            u = gpcov(pivot[i], pivot[k], covpar, include);
            for(R[k][i] = u, l = 0; l < k; l++)
                R[k][i] -= R[l][i] * R[l][k];
            R[k][i] /= R[k][k];
            d[i] -= R[k][i] * R[k][i];
        }
        
        k++;
        flag = (k < max_rank);
        if(flag && dostopping){
            for(max_diag = k, i = k + 1; i < max_scan; i++)
                if(d[i] > d[max_diag])
                    max_diag = i;
            flag = (d[max_diag] > tol);
        }
    }
    
    rank[0] = k;
    if(padzero){
        for(l = k; l < N; l++)
            d[l] = 0.0;
    }
}

void spchol_agp(double **R, int N, double tol, int *pivot, int *rank, int max_rank, double gpcov(int, int, double*, int*), int ncomp, int npar, int nvar, double *covpar, int *include, int max_scan, double *d, int dopivoting, int dostopping, int padzero, double *lpen, int *active){
    
    
    // Sparse cholesky factorization with pivoting and diagonal augmentation
    // for stationary covariance functions.
    // Accepts an empty matrix R and a function gpcov(i, j, ...) that is called
    // to compute the (i,j)-th element of the original covariance matrix
    
    set_lower_tri_zero(R, N, max_rank);
    
    int i, a, l, shift1, shift2;
    double u, b;
    int randFirst = 1;
    for(i = 0; i < N; i++) pivot[i] = i;
    if(dopivoting) q_sort(lpen, pivot, max_scan, randFirst);
    
    for(i = 0; i < N; i++) {
        d[i] = 0.0; shift1 = 0; shift2 = 0;
        for(l = 0; l < ncomp; l++){
            if(active[l]) d[i] += gpcov(pivot[i], pivot[i], covpar + shift1, include + shift2);
            shift1 += npar;
            shift2 += nvar;
        }
    }
    
    int k = 0, max_diag;
    for(max_diag = k, i = k + 1; i < max_scan; i++)
        if(d[i] > d[max_diag])
            max_diag = i;
    
    tol *= d[max_diag];
    int flag = (k < max_rank);
    if(dostopping)
        flag = (d[max_diag] > tol);
    
    while(flag){
        if(dopivoting){
            if(max_diag > k){
                a = pivot[k];
                pivot[k] = pivot[max_diag];
                pivot[max_diag] = a;
                
                b = d[k];
                d[k] = d[max_diag];
                d[max_diag] = b;
                
                for(i = 0; i < k; i++){
                    b = R[i][k];
                    R[i][k] = R[i][max_diag];
                    R[i][max_diag] = b;
                }
            }
        }
        
        R[k][k] = sqrt(d[k]);
        
        for(i = k + 1; i < N; i++){
            u = 0.0; shift1 = 0; shift2 = 0;
            for(l = 0; l < ncomp; l++){
                if(active[l]) u += gpcov(pivot[i], pivot[k], covpar + shift1, include + shift2);
                shift1 += npar;
                shift2 += nvar;
            }
            for(R[k][i] = u, l = 0; l < k; l++)
                R[k][i] -= R[l][i] * R[l][k];
            R[k][i] /= R[k][k];
            d[i] -= R[k][i] * R[k][i];
        }
        
        k++;
        flag = (k < max_rank);
        if(flag && dostopping){
            for(max_diag = k, i = k + 1; i < max_scan; i++)
                if(d[i] > d[max_diag])
                    max_diag = i;
            flag = (d[max_diag] > tol);
        }
    }
    
    rank[0] = k;
    if(padzero){
        for(l = k; l < N; l++)
            d[l] = 0.0;
    }
}

void chol(double **R, int N, double tol, int *pivot, int *rank, int max_rank,
          double *d, double **A, int dopivoting, int padzero){
    
    set_lower_tri_zero(R, N, max_rank);
    
    int i, a, l;
    double u, b;
    
    for(i = 0; i < N; i++){
        pivot[i] = i;
        d[i] = A[i][i];
    }
    
    int k = 0, max_diag;
    for(max_diag = k, i = k + 1; i < N; i++)
        if(d[i] > d[max_diag])
            max_diag = i;
    int flag = (d[max_diag] > tol);
    
    
    while(flag){
        if(dopivoting){
            if(max_diag > k){
                a = pivot[k];
                pivot[k] = pivot[max_diag];
                pivot[max_diag] = a;
                
                b = d[k];
                d[k] = d[max_diag];
                d[max_diag] = b;
                
                for(i = 0; i < k; i++){
                    b = R[i][k];
                    R[i][k] = R[i][max_diag];
                    R[i][max_diag] = b;
                }
            }
        }
        
        R[k][k] = sqrt(d[k]);
        
        for(i = k + 1; i < N; i++){
            u = A[pivot[i]][pivot[k]];
            for(R[k][i] = u, l = 0; l < k; l++)
                R[k][i] -= R[l][i] * R[l][k];
            R[k][i] /= R[k][k];
            d[i] -= R[k][i] * R[k][i];
        }
        
        k++;
        flag = (k < max_rank);
        if(flag){
            for(max_diag = k, i = k + 1; i < N; i++)
                if(d[i] > d[max_diag])
                    max_diag = i;
            flag = (d[max_diag] > tol);
        }
    }
    
    rank[0] = k;
    if(padzero){
        for(l = k; l < N; l++)
            d[l] = 0.0;
    }
}

void trisolve(double **R, int m, double *b, double *x, int transpose){
    
    int i, j;
    if(transpose){
        for(j = 0; j < m; j++){
            for(x[j] = b[j], i = 0; i < j; i++)
                x[j] -= x[i] * R[i][j];
            x[j] /= R[j][j];
        }
    } else {
        for(j = m - 1; j >= 0; j--){
            for(x[j] = b[j], i = j + 1; i < m; i++)
                x[j] -= R[j][i] * x[i];
            x[j] /= R[j][j];
        }
    }
}

void triprod(double **R, int m, int n, double *x, double *b, int transpose){
    
    int i, j;
    if(transpose){
        for(i = 0; i < m; i++)
            for(b[i] = 0.0, j = 0; j <= i; j++)
                b[i] += R[j][i] * x[j];
        for(; i < n; i++)
            for(b[i] = 0.0, j = 0; j < m; j++)
                b[i] += R[j][i] * x[j];
    } else{
        for(i = 0; i < m; i++)
            for(b[i] = 0.0, j = i; j < n; j++)
                b[i] += R[i][j] * x[j];
    }
}

// GLOBAL CONSTANTS
const double rho_min_trunc = 0.0000000001;

// GLOBAL VARIABLES
int n, nnew, p, dmax, max_rank, toprint, dopivot, nlam, nrho, ngrid;
double **x, *y, **lpenMat, *logpmodel, sig_shape, sig_rate, *lsgrid, *lamsq, *rhosq, *lplam, *lprho;

// VARIABLE REPEATEDLY USED BY AUX FUNCTIONS
int *rank, *pivot, *grank, *gpvt, *activestatus;
double *d, *a, *gd, *yp, *z, *ryp, *gptol, *cpar, *lpen, *v, *u, *mu_hat, *znorm, *Hz, *Rz;
double **R, **F, **G, *Rad, **GF, **H, **HC;


// Main auxiliary functions
double sevscov(int i, int j, double *covpar, int *include);
void lsFn(int *include, double *yRes, int ilam, double *ls_elements, int maximum_rank, double *sig_pars);
double works(int *include, int csize, double *yRes, double *par, int maximum_rank, double pactive, double *sig_pars);
double rgpFn(int *include, int csize, double *yRes, double *covpar, double *fdraw, double *sig_pars);
double ragpFn(int ncomp, int npar, int *include, double *yRes, double *covpar, double *fmean, double *fdraw);

// Codes for sparse additive-interactive regression with GP




double sevscov(int i, int j, double *covpar, int *include){
    
    int l, nact = 0;
    double dist = 0.0, diff;
    for(l = 0; l < p; l++){
        if(include[l]){
            nact++;
            diff = x[i][l] - x[j][l];
            dist += diff * diff;
        }
    }
    if(nact == 0) nact = 1;
    //nact = 1;
    double bsq = exp(log(covpar[1])/(double)nact);
    return covpar[0] * exp(-bsq * dist);
}

// Log-likelihood score calculator. Depends on the following
// (A) global variables: *cpar, **R, n, *gptol, *pivot, *rank, *d, *a, dopivot, *lpen, **Rad, **F, **G, *gpvt, *grank, *gd, *yp, *ryp, *z, *lamsq, *rhosq,, *lplam, *lprho
// (B) functions: spchol_st, triprod, chol, trisolve, sumsquares,

void lsFn(int *include, double *yRes, int ilam, double *ls_elements, int maximum_rank, double *sig_pars){
    
    int i, irho, j;
    double SSR = 0.0, logD = 0.0;
    cpar[0] = 1.0; cpar[1] = lamsq[ilam];
    spchol_st(R, n, gptol[0], pivot, rank, maximum_rank, sevscov, cpar, include, n, d, dopivot, 1, 0, lpen);
    
    for(irho = 0; irho < nrho; irho++){
        if(rhosq[irho] > rho_min_trunc){
            for(i = 0; i < rank[0]; i++) a[pivot[i]] = 1.0 / rhosq[irho];
            for(i = rank[0]; i < n; i++) a[pivot[i]] = 1.0 / rhosq[irho] + d[i];
            
            for(i = 0; i < rank[0]; i++){
                for(j = 0; j < n; j++)
                    Rad[j] = R[i][j] / a[pivot[j]];
                triprod(R, rank[0], n, Rad, F[i], 0);
            }
            for(i = 0; i < rank[0]; i++) F[i][i] += 1.0;
            
            chol(G, rank[0], gptol[1], gpvt, grank, rank[0], gd, F, 0, 0);
            
            double logdetA = 0.0, logdetG = 0.0;
            
            for(i = 0; i < n; i++) logdetA += log(a[i]);
            for(i = 0; i < rank[0]; i++) logdetG += log(G[i][i]);
            
            for(i = 0; i < n; i++) yp[i] = yRes[pivot[i]] / a[pivot[i]];
            triprod(R, rank[0], n, yp, ryp, 0);
            trisolve(G, rank[0], ryp, z, 1);
            
            for(i = 0; i < n; i++) yp[i] = yRes[i] / sqrt(a[i]);
            SSR = (sumsquares(yp,n) - sumsquares(z,rank[0])) / rhosq[irho];
            logD = logdetA + 2.0 * logdetG + (double)n * log(rhosq[irho]);
            
        } else {
            SSR = sumsquares(yRes, n);
            logD = 0.0;
        }
        ls_elements[irho] = -0.5*logD - (0.5*(double)n + sig_pars[0]) * log(0.5*SSR + sig_pars[1]) + lplam[ilam] + lprho[irho];
    }
}

double works(int *include, int csize, double *yRes, double *par, int maximum_rank, double pactive, double *sig_pars){
    int i, ilam, j, skip = 0;
    double ls = 0.0;
    
    if(csize > 0){
        vect_set(lpen, n, 0.0);
        for(j = 0; j < p; j++) if(include[j]) for(i = 0; i < n; i++) lpen[i] += lpenMat[j][i];
        for(ilam = 0; ilam < nlam; ilam++){
            lsFn(include, yRes, ilam, lsgrid + skip, maximum_rank, sig_pars);
            skip += nrho;
        }
        ls = logsum(lsgrid, ngrid) + logpmodel[csize-1] + log(pactive);
        int ii = rdraw(ngrid, lsgrid, 1);
        par[0] = rhosq[ii % nrho]; par[1] = lamsq[ii / nrho];
    } else {
        ls = - (0.5*(double)n + sig_pars[0]) * log(0.5*sumsquares(yRes,n) + sig_pars[1]) + log(1.0 - pactive);
        par[0] = 0.0; par[1] = 0.0;
    }
    return ls;
}

double rgpFn(int *include, int csize, double *yRes, double *covpar, double *fdraw, double *sig_pars){
    
    int i, j;
    double sigma, SSR;

    if(csize == 0 || covpar[0] == 0.0){
        SSR = sumsquares(yRes, n);
        sigma = sqrt((0.5*SSR + sig_pars[1]) / rgamma(0.5*(double)n + sig_pars[0], 1.0));
        for(i = 0; i < n; i++) fdraw[i] = 0.0;
    } else {
        spchol_st(R, n, gptol[0], pivot, rank, max_rank, sevscov, covpar, include, n, d, dopivot, 1, 0, lpen);
        
        for(i = 0; i < rank[0]; i++) a[pivot[i]] = 1.0;
        for(i = rank[0]; i < n; i++) a[pivot[i]] = 1.0 + d[i];
        
        for(i = 0; i < rank[0]; i++){
            for(j = 0; j < n; j++)
                Rad[j] = R[i][j] / a[pivot[j]];
            triprod(R, rank[0], n, Rad, F[i], 0);
        }
        for(i = 0; i < rank[0]; i++) F[i][i] += 1.0;
        chol(G, rank[0], 0.0, gpvt, grank, rank[0], gd, F, 0, 0);
        for(i = 0; i < rank[0]; i++) F[i][i] -= 1.0;		// back to F0
        
        for(i = 0; i < n; i++) yp[i] = yRes[pivot[i]] / a[pivot[i]];
        triprod(R, rank[0], n, yp, ryp, 0);
        trisolve(G, rank[0], ryp, z, 1);
        trisolve(G, rank[0], z, v, 0);
        for(i = 0; i < rank[0]; i++) u[i] = ryp[i] - inprod(F[i], v, rank[0]);
        
        triprod(R, rank[0], n, u, mu_hat, 1);
        
        for(i = 0; i < rank[0]; i++) trisolve(G, rank[0], F[i], GF[i], 1);
        for(i = 0; i < rank[0]; i++){
            for(j = 0; j < rank[0]; j++) H[i][j] = inprod(GF[i], GF[j], rank[0]) - F[i][j];
            H[i][i] += 1.0;
        }
        chol(HC, rank[0], 0.0, gpvt, grank, rank[0], gd, H, 0, 0);
        for(i = 0; i < rank[0]; i++) znorm[i] = rnorm(0.0, 1.0);
        triprod(HC, rank[0], rank[0], znorm, Hz, 1);
        triprod(R, rank[0], n, Hz, Rz, 1);
        
        for(i = 0; i < n; i++) yp[i] = yRes[i] / sqrt(a[i]);
        SSR = sumsquares(yp,n) - sumsquares(z,rank[0]);
        sigma = sqrt((0.5*SSR + sig_pars[1]) / rgamma(0.5*(double)n + sig_pars[0], 1.0));
        for(i = 0; i < n; i++) fdraw[pivot[i]] = mu_hat[i] + sigma * Rz[i];
    }
    return sigma;
}

// Functions needeed for PREDICTION //
double ragpFn(int ncomp, int npar, int *include, double *yRes, double *covpar, double *fmean, double *fdraw){
    
    int i, j, l, N = n + nnew;
    double sigma = 0.0, totsignal = 0.0, SSR = 0.0;
    
    for(l = 0; l < ncomp; l++){
        activestatus[l] = any(include + l*p, p);
        if(activestatus[l]) totsignal += covpar[npar*l];
    }
    
    if(totsignal == 0.0){
        SSR = sumsquares(yRes, n);
        sigma = sqrt((0.5*SSR + sig_rate) / rgamma(0.5*(double)n + sig_shape, 1.0));
        for(i = 0; i < N; i++) fmean[i] = fdraw[i] = 0.0;
    } else {
        for(i = 0; i < n; i++){
            lpen[i] = 0.0;
            for(l = 0; l < ncomp; l++) if(activestatus[l]) for(j = 0; j < p; j++) if(include[l*p + j]) lpen[i] += sqrt(covpar[npar*l])*lpenMat[j][i];
        }
        spchol_agp(R, N, gptol[0], pivot, rank, max_rank, sevscov, ncomp, npar, p, covpar, include, n, d, dopivot, 1, 0, lpen, activestatus);
        
        
        for(i = 0; i < rank[0]; i++) a[pivot[i]] = 1.0;
        for(i = rank[0]; i < n; i++) a[pivot[i]] = 1.0 + d[i];
        
        for(i = 0; i < rank[0]; i++){
            for(j = 0; j < n; j++)
                Rad[j] = R[i][j] / a[pivot[j]];
            triprod(R, rank[0], n, Rad, F[i], 0);
        }
        for(i = 0; i < rank[0]; i++) F[i][i] += 1.0;
        chol(G, rank[0], 0.0, gpvt, grank, rank[0], gd, F, 0, 0);
        for(i = 0; i < rank[0]; i++) F[i][i] -= 1.0;		// back to F0
        
        for(i = 0; i < n; i++) yp[i] = yRes[pivot[i]] / a[pivot[i]];
        triprod(R, rank[0], n, yp, ryp, 0);
        trisolve(G, rank[0], ryp, z, 1);
        trisolve(G, rank[0], z, v, 0);
        for(i = 0; i < rank[0]; i++) u[i] = ryp[i] - inprod(F[i], v, rank[0]);
        
        triprod(R, rank[0], N, u, mu_hat, 1);
        
        for(i = 0; i < rank[0]; i++) trisolve(G, rank[0], F[i], GF[i], 1);
        for(i = 0; i < rank[0]; i++){
            for(j = 0; j < rank[0]; j++) H[i][j] = inprod(GF[i], GF[j], rank[0]) - F[i][j];
            H[i][i] += 1.0;
        }
        chol(HC, rank[0], 0.0, gpvt, grank, rank[0], gd, H, 0, 0);
        for(i = 0; i < rank[0]; i++) znorm[i] = rnorm(0.0, 1.0);
        triprod(HC, rank[0], rank[0], znorm, Hz, 1);
        triprod(R, rank[0], N, Hz, Rz, 1);
        
        for(i = 0; i < n; i++) yp[i] = yRes[i] / sqrt(a[i]);
        SSR = sumsquares(yp,n) - sumsquares(z,rank[0]);
        sigma = sqrt((0.5*SSR + sig_rate) / rgamma(0.5*(double)n + sig_shape, 1.0));
        
        for(i = 0; i < N; i++){
            fmean[pivot[i]] = mu_hat[i];
            fdraw[pivot[i]] = mu_hat[i] + sigma * Rz[i];
        }
    }
    return sigma;
}


// Interface functions *** not declared ***

void airGP_one(double *xvar, double *yvar, int *dim, double *lamsqR, double *rhosqR, double *hpar, double *logp_modelsize, double *lplamR, double *lprhoR, double *lpenalty, double *ls1, double *ls0, double *tolchol){
    
    n = dim[0];
    p = dim[1];
    max_rank = dim[2];
    nlam = dim[3];
    nrho = dim[4];
    toprint = 1;
    dopivot = 1;
    
    ngrid = nlam * nrho;
    lamsq = lamsqR;
    rhosq = rhosqR;
    lplam = lplamR;
    lprho = lprhoR;
    
    sig_shape = hpar[0];
    sig_rate = hpar[1];
    logpmodel = logp_modelsize;
    //Rprintf("logpmodel[0] = %g\n", logpmodel[0]);
    
    int i, j;
    
    x = mymatrix(n, p);
    //y = yvar;
    R = mymatrix(max_rank, n);
    d = vect(n);
    a = vect(n);
    rank = ivect(1); rank[0] = max_rank;
    pivot = ivect(n);
    F = mymatrix(max_rank, max_rank);
    G = mymatrix(max_rank, max_rank);
    Rad = vect(n);
    gd = vect(max_rank);
    grank = ivect(1);
    gpvt = ivect(max_rank);
    yp = vect(n);
    z = vect(max_rank);
    ryp = vect(max_rank);
    gptol = tolchol;
    cpar = vect(2);
    lpenMat = mymatrix(p, n);
    lpen = vect(n);
    lsgrid = vect(ngrid);
    
    v = vect(max_rank); u = vect(max_rank); mu_hat = vect(n); znorm = vect(max_rank); Hz = vect(max_rank), Rz = vect(n);
    GF = mymatrix(max_rank, max_rank); H = mymatrix(max_rank, max_rank); HC = mymatrix(max_rank, max_rank);

    int pos = 0;
    for(j = 0; j < p; j++){
        for(i = 0; i < n; i++) x[i][j] = xvar[pos  + i];
        pos += n;
    }
    
    pos = 0;
    for(j = 0; j < p; j++){
        for(i = 0; i < n; i++) lpenMat[j][i] = lpenalty[pos + i];
        pos += n;
    }
    
    GetRNGstate();
    double *par_one = vect(2); vect_set(par_one, 2, 0.0);
    int *ix = ivect(p);
    for(j = 0; j < p; j++) ix[j] = 0;
    ls0[0] = works(ix, 0, yvar, par_one, max_rank, 0.5, hpar);
    //Rprintf("ls0 = %g\n", ls0[0]);
    for(j = 0; j < p; j++){
        ix[j] = 1;
        ls1[j] = works(ix, 1, yvar, par_one, max_rank, 0.5, hpar);
        ix[j] = 0;
    }
    PutRNGstate();
    //Rprintvec("ls1 =", "%g ", ls1, p);
    
}


// MAIN FUNCTION FOR FITTING THE ADDITIVE GP MODEL by the TOGGLE METHOD
void airGP_tog(double *data, int *state, double *state_pars, double *fcompvec, int *comp_size, char **compix, int *dim, int *nhbr_ix, int *nhbr_len, double *nhbr_val, double *gpgrid, double *hpar, double *comp_adj, double *logp_modelsize, double *pmove, double *move_diagnostic, double *lpenalty, double *varp, char **ix_store, double *par_store, double *sigma_store, double *ftot_store, double *pactive_store, double *rhob_store, double *tolchol, int *print_details){
    
    //Rprintf("Inside C code\n");
    
    int i, j, k, irho;
    //Rprintf("rho_min_trunc = %g\n", rho_min_trunc);
    double rho_min_trunc_2 = 0.0;

    n = dim[0];
    p = dim[1];
    int ncomp = dim[2];
    max_rank = dim[3];
    nlam = dim[4];
    nrho = dim[5];
    dmax = dim[6];
    int nsamp = dim[7];
    int thin = dim[8];
    //int nsweep = nsamp*thin;
    
    toprint = 1;
    dopivot = 1;
    
    ngrid = nlam * nrho;
    lamsq = vect(nlam); vect_import(lamsq, nlam, gpgrid);
    lplam = vect(nlam); vect_import(lplam, nlam, gpgrid + nlam);
    rhosq = vect(nrho); vect_set(rhosq, nrho, 0.0);
    lprho = vect(nrho); vect_set(lprho, nrho, -log((double)nrho));
    
    sig_shape = hpar[0];
    sig_rate = hpar[1];
    double pactive_a = hpar[2], pactive_b = hpar[3], rho_a = hpar[4], rho_b_shape = hpar[5], rho_b_Hmean = hpar[6];
    logpmodel = logp_modelsize;

    x = mymatrix(n, p);
    y = data;
    R = mymatrix(max_rank, n);
    d = vect(n);
    a = vect(n);
    rank = ivect(1); rank[0] = max_rank;
    pivot = ivect(n);
    F = mymatrix(max_rank, max_rank);
    G = mymatrix(max_rank, max_rank);
    Rad = vect(n);
    gd = vect(max_rank);
    grank = ivect(1);
    gpvt = ivect(max_rank);
    yp = vect(n);
    z = vect(max_rank);
    ryp = vect(max_rank);
    gptol = tolchol;
    cpar = vect(2);
    lpenMat = mymatrix(p, n);
    lpen = vect(n);
    lsgrid = vect(ngrid);
    activestatus = ivect(ncomp);
    
    v = vect(max_rank); u = vect(max_rank); mu_hat = vect(n); znorm = vect(max_rank); Hz = vect(max_rank), Rz = vect(n);
    GF = mymatrix(max_rank, max_rank); H = mymatrix(max_rank, max_rank); HC = mymatrix(max_rank, max_rank);
    
    double *xvar = data + n;
    int pos = 0;
    for(j = 0; j < p; j++) {
        for(i = 0; i < n; i++) x[i][j] = xvar[pos  + i]; pos += n;
    }
    
    pos = 0;
    for(j = 0; j < p; j++){
        for(i = 0; i < n; i++) lpenMat[j][i] = lpenalty[pos + i]; pos += n;
    }
    
    double pactive = state_pars[0];
    double sigma = state_pars[1];
    double sigmasq = sigma*sigma;
    double rho_b = state_pars[2];
    double *par = state_pars + 3;
    
    int *ixlong = ivect(p*ncomp);
    int **ix = (int **)R_alloc(ncomp, sizeof(int *));
    for(k = 0; k < ncomp; k++) ix[k] = ixlong + k*p;
    
    double **fcomp = (double **)R_alloc(ncomp, sizeof(double *));
    for(k = 0; k < ncomp; k++) fcomp[k] = fcompvec + k*n;
    double *tausq_comp = vect(ncomp);
    vect_set(tausq_comp, ncomp, 0.0);
    
    char b[10000];
    for(k = 0; k < ncomp; k++){
        for(j = 0; j < p; j++) ix[k][j] = 0;
        strcpy(b, compix[k]);
        string_locator(b, ix[k], p);
    }
    double *pp = vect(p), ls;

    GetRNGstate();
    double *yRes = vect(n), *ftot = vect(n), *sig_pars = vect(2);
    vect_import(yRes, n, y); vect_set(ftot, n, 0.0);
    double *rhoAfactor, *rhoBfactor;
    rhoAfactor = comp_adj; rhoBfactor = comp_adj + ncomp;
    
    if(state[0]){
        //Rprintf("Starting state initialization...");
        sig_pars[0] = sig_shape; sig_pars[1] = sig_rate;
        double base = 0.00001 * vmax(varp, p);
        if(base < 0.00001) base = 0.00001;
        for(j = 0; j < p; j++) pp[j] = base + varp[j];
        for(i = 0; i < dmax; i++) logpmodel[i] += lchoose(p, 1+i);
        for(k = 0; k < ncomp; k++){
            comp_size[k] = 0; activestatus[k] = 0; tausq_comp[k] = 0.0;
            for(j = 0; j < p; j++) ix[k][j] = 0;
            activestatus[k] = (runif(0.0, 1.0) < pactive);
            if(activestatus[k]){
                comp_size[k] = 1 + rdraw(dmax, logpmodel, 1);
                for(i = 0; i < comp_size[k]; i++){
                    j = rdraw(p, pp, 0);
                    ix[k][j] = 1;
                    pp[j] = base;
                }
                for(irho = 0; irho < nrho; irho++) rhosq[irho] = stabler(rgamma(rho_a*rhoAfactor[k], 1.0)/(rho_b*rhoBfactor[k]), rho_min_trunc_2);
            }
            ls = works(ix[k], comp_size[k], yRes, par + 2*k, max_rank, pactive, sig_pars);
            sigma = rgpFn(ix[k], comp_size[k], yRes, par + 2*k, fcomp[k], sig_pars);
            sigmasq = sigma*sigma;
            if(activestatus[k]){
                tausq_comp[k] = sigmasq * par[2*k];
                for(i = 0; i < n; i++){
                    ftot[i] += fcomp[k][i];
                    yRes[i] -= fcomp[k][i];
                }
                sig_pars[0] += rho_a * rhoAfactor[k];
                sig_pars[1] += tausq_comp[k] * rho_b * rhoBfactor[k];
            }
        }
        for(i = 0; i < dmax; i++) logpmodel[i] -= lchoose(p, 1+i);
        //Rprintf("Done\n");
    } else {
        for(k = 0; k < ncomp; k++){
            activestatus[k] = comp_size[k] ? 1 : 0;
            if(activestatus[k]){
                tausq_comp[k] = sigmasq * par[2*k];
                for(i = 0; i < n; i++) ftot[i] += fcomp[k][i];
            }
        }
    }
    //Rprintvec("par = ", "%g ", par, 2*ncomp);
    
    
    int nActive = isum(activestatus, ncomp);
    int numcall_works = 0;
    int samp, reps, sweep, move_type, ticker = (int)ceil((double)nsamp / 10.0);
    double *nprop = move_diagnostic, *nacpt = move_diagnostic + 4;
    nprop[0] = nprop[1] = nprop[2] = nprop[3] = 0.0;
    nacpt[0] = nacpt[1] = nacpt[2] = nacpt[3] = 0.0;
    double *varpc = vect(p), *varp_k = vect(p), *varpc_k = vect(p);
    vect_set(varpc, p, 1.0);
    
    int jj, jj2, f_pos = 0, par_pos = 0, ix_pos = 0, nhbr_skip, *ixnew = ivect(p);
    double uu, log_alpha, pfwd, pbwd, lsnew, incr, *parnew = vect(2);
    double pactive_acpt_logprob, pactive_prop;
    int pactive_nacpt = 0;
    
    double *rtimes = vect(3), Rsqk = 0.0;
    time_t t0, t1;
    rtimes[0] = 0.0; rtimes[1] = 0.0; rtimes[2] = 0.0;
    pmove[1] = 0.0; pmove[2] = 0.0;
    
    double temp = 1.0, SSR = 0.0;
    if(print_details[0]) Rprintf("airGP: n = %d, p = %d, ncomp = %d, E(nActive) = %.1f\n", n, p, ncomp, (double)ncomp*(pactive));
    
    sig_pars[0] = sig_shape + rho_a * sum_subset(rhoAfactor, ncomp, activestatus);
    sig_pars[1] = sig_rate + rho_b * inprod_subset(tausq_comp, rhoBfactor, ncomp, activestatus);
    
    sweep = 0;
    for(samp = 0; samp < nsamp; samp++){
        for(reps = 0; reps < thin; reps++){
            sweep++;
            
            for(k = 0; k < ncomp; k++){
                for(irho = 0; irho < nrho; irho++) rhosq[irho] = stabler(rgamma(rho_a * rhoAfactor[k], 1.0)/(rho_b * rhoBfactor[k]), rho_min_trunc_2);
                if(activestatus[k]){
                    for(i = 0; i < n; i++) ftot[i] -= fcomp[k][i];
                    sig_pars[0] -= rho_a * rhoAfactor[k];
                    sig_pars[1] -= tausq_comp[k] * rho_b * rhoBfactor[k];
                    rhosq[0] = stabler(tausq_comp[k] / sigmasq, rho_min_trunc_2);
                }

                for(i = 0; i < n; i++) yRes[i] = y[i] - ftot[i];
                ls = works(ix[k], comp_size[k], yRes, par + 2*k, max_rank, pactive, sig_pars);
                numcall_works++;
                
                move_type = rdraw(4, pmove + 4 * comp_size[k], 0);
                nprop[move_type] += 1.0;
                
                t0 = clock();
                switch(move_type){
                        
                    case 0 : // add
                        for(j = 0; j < p; j++){
                            varp_k[j] = ix[k][j] ? 0.0 : varp[j];
                            ixnew[j] = ix[k][j];
                        }
                        jj = rdraw(p, varp_k, 0);
                        ixnew[jj] = 1;
                        lsnew = works(ixnew, comp_size[k]+1, yRes, parnew, max_rank, pactive, sig_pars);
                        numcall_works++;
                        pfwd = pmove[4 * comp_size[k]] * varp_k[jj] / sum(varp_k, p);
                        
                        for(j = 0; j < p; j++) varpc_k[j] = ixnew[j] ? varpc[j] : 0.0;
                        pbwd = pmove[4 * comp_size[k] + 4 + 1] * varpc_k[jj] / sum(varpc_k, p);
                        
                        log_alpha = lsnew - ls + log(pbwd) - log(pfwd);
                        uu = runif(0.0, 1.0);
                        if(log(uu) < temp * log_alpha){
                            nacpt[0] += 1.0;
                            ix[k][jj] = 1;
                            par[2*k] = parnew[0]; par[2*k + 1] = parnew[1];
                            comp_size[k]++;
                        }
                        break;
                        
                    case 1 : //remove
                        for(j = 0; j < p; j++){
                            varpc_k[j] = ix[k][j] ? varpc[j] : 0.0;
                            ixnew[j] = ix[k][j];
                        }
                        jj = rdraw(p, varpc_k, 0);
                        ixnew[jj] = 0;
                        lsnew = works(ixnew, comp_size[k]-1, yRes, parnew, max_rank, pactive, sig_pars);
                        numcall_works++;
                        pfwd = pmove[4 * comp_size[k] + 1] * varpc_k[jj] / sum(varpc_k, p);
                        
                        for(j = 0; j < p; j++) varp_k[j] = ixnew[j] ?  0.0 : varp[j];
                        pbwd = pmove[4 * comp_size[k] - 4] * varp_k[jj] / sum(varp_k, p);
                        
                        log_alpha = lsnew - ls + log(pbwd) - log(pfwd);
                        uu = runif(0.0, 1.0);
                        if(log(uu) < temp * log_alpha){
                            nacpt[1] += 1.0;
                            ix[k][jj] = 0;
                            par[2*k] = parnew[0]; par[2*k + 1] = parnew[1];
                            comp_size[k]--;
                        }
                        break;
                        
                    case 2 : //swap
                        for(j = 0; j < p; j++){
                            varp_k[j] = ix[k][j] ? 0.0 : varp[j];
                            varpc_k[j] = ix[k][j] ? varpc[j] : 0.0;
                            ixnew[j] = ix[k][j];
                        }
                        jj = rdraw(p, varp_k, 0);
                        jj2 = rdraw(p, varpc_k, 0);
                        ixnew[jj] = 1; ixnew[jj2] = 0;
                        lsnew = works(ixnew, comp_size[k], yRes, parnew, max_rank, pactive, sig_pars);
                        numcall_works++;
                        pfwd = (varp_k[jj] / sum(varp_k, p)) * (varpc_k[jj2] / sum(varpc_k, p));
                        varp_k[jj] = 0.0; varp_k[jj2] = varp[jj2];
                        varpc_k[jj] = varpc[jj]; varpc_k[jj2] = 0.0;
                        pbwd = (varpc_k[jj] / sum(varpc_k, p)) * (varp_k[jj2] / sum(varp_k, p));
                        
                        log_alpha = lsnew - ls + log(pbwd) - log(pfwd);
                        uu = runif(0.0, 1.0);
                        if(log(uu) < temp * log_alpha){
                            nacpt[2] += 1.0;
                            ix[k][jj] = 1; ix[k][jj2] = 0;
                            par[2*k] = parnew[0]; par[2*k + 1] = parnew[1];
                        }
                        break;
                        
                    default : //refresh
                        nacpt[3] += 1.0;
                        break;
                }
                t1 = clock();
                rtimes[0] += (double)(t1 - t0) / CLOCKS_PER_SEC;

                t0 = clock();
                sigma = rgpFn(ix[k], comp_size[k], yRes, par + 2*k, fcomp[k], sig_pars);
                sigmasq = sigma * sigma;

                tausq_comp[k] = 0.0;
                activestatus[k] = comp_size[k] ? 1 : 0;
                if(activestatus[k]){
                    tausq_comp[k] = par[2*k] * sigmasq;
                    for(i = 0; i < n; i++) ftot[i] += fcomp[k][i];
                    //sig_pars[0] += rho_a/comp_adj;
                    sig_pars[0] += rho_a * rhoAfactor[k];
                    sig_pars[1] += tausq_comp[k] * rho_b * rhoBfactor[k];
                }

                t1 = clock();
                rtimes[1] += (double)(t1 - t0) / CLOCKS_PER_SEC;
                
                
                t0 = clock();
                incr = R_pow((double)(1+sweep), -0.67);
                if(incr > 0.1) incr = 0.1;
                
                if(activestatus[k]){
                    Rsqk = sigmoid_lin(tausq_comp[k]/sigmasq);
                    nhbr_skip = 0;
                    for(jj = 0; jj < p; jj++){
                        if(ix[k][jj]) {
                            for(j = 0; j < nhbr_len[jj]; j++) varp[nhbr_ix[nhbr_skip + j]] += Rsqk * incr * nhbr_val[nhbr_skip + j];
                        }
                        nhbr_skip += nhbr_len[jj];
                    }
                }
                //for(j = 0; j < p; j++) varpc[j] = 1.0; // / varp[j];
                
                t1 = clock();
                rtimes[2] += (double)(t1 - t0) / CLOCKS_PER_SEC;
                
            }

            for(i = 0; i < n; i++) yRes[i] = y[i] - ftot[i];
            SSR = sumsquares(yRes, n);
            sigma = sqrt((0.5*SSR + sig_pars[1]) / rgamma(0.5*(double)n + sig_pars[0], 1.0));
            sigmasq = sigma * sigma;

            nActive = isum(activestatus, ncomp);
            pactive_prop = rbeta(pactive_a + (double)nActive, pactive_b + (double)(ncomp - nActive));
            pactive_acpt_logprob = 0.0;
            if(log(runif(0.0, 1.0)) < pactive_acpt_logprob) {
                pactive_nacpt++;
                pactive = pactive_prop;
            }
            
            rho_b = rgamma(rho_a * sum_subset(rhoAfactor, ncomp, activestatus) + rho_b_shape, 1.0)/(inprod_subset(tausq_comp, rhoBfactor, ncomp, activestatus) / sigmasq + (rho_b_shape - 1.0)/rho_b_Hmean);

            sig_pars[0] = sig_shape + rho_a * sum_subset(rhoAfactor, ncomp, activestatus);
            sig_pars[1] = sig_rate + rho_b * inprod_subset(tausq_comp, rhoBfactor, ncomp, activestatus);
        }

        // Storage and printing
        pactive_store[samp] = pactive;
        sigma_store[samp] = sigma;
        rhob_store[samp] = rho_b;
        for(k = 0; k < ncomp; k++){
            par_store[par_pos] = tausq_comp[k] / sigmasq; par_pos++;
            par_store[par_pos] = par[2*k + 1]; par_pos++;
            locator_string(ix[k], p, ix_store[ix_pos]); ix_pos++;
        }
        for(i = 0; i < n; i++) ftot_store[f_pos + i] = ftot[i]; f_pos += n;
        
        if(print_details[0]){
            if((samp + 1) % ticker == 0){
                Rprintf("samp = %d. Active: ", samp + 1);
                for(k = 0; k < ncomp; k++) if(comp_size[k] > 0) Rprintf(" %s ", ix_store[ix_pos - ncomp + k]);
                Rprintf(" nActive: %d\n", nActive);
            }
        }
    }
    
    PutRNGstate();
    //Rprintf("pactive acceptance rate = %g\n", (double)pactive_nacpt / (double)sweep);
    state_pars[0] = pactive;
    state_pars[1] = sigma;
    state_pars[2] = rho_b;
    for(k = 0; k < ncomp; k++) locator_string(ix[k], p, compix[k]);

}

//PREDICT function
void airGP_pred(double *xvar, double *yvar, double *xnew, int *dim, double *hpar, double *lpenalty, char **ix_store, double *par_store, double *tolchol, double *fmeanstore, double *fstore, double *sigstore){
    
    n = dim[0];
    p = dim[1];
    int ncomp = dim[2];
    max_rank = dim[3];
    nnew = dim[4];
    int nsamp = dim[5];
    toprint = 1;
    dopivot = 1;
    
    sig_shape = hpar[0];
    sig_rate = hpar[1];
    
    //Rprintf("n = %d, p = %d, nnew = %d, ncomp = %d, nsamp = %d\n", n, p, nnew, ncomp, nsamp);
    int i, j, l, N = n + nnew;
    
    x = mymatrix(n + nnew, p);
    y = yvar;
    R = mymatrix(max_rank, n + nnew);
    d = vect(n + nnew);
    a = vect(n + nnew);
    rank = ivect(1); rank[0] = max_rank;
    pivot = ivect(n + nnew);
    F = mymatrix(max_rank, max_rank);
    G = mymatrix(max_rank, max_rank);
    Rad = vect(n);
    gd = vect(max_rank);
    grank = ivect(1);
    gpvt = ivect(max_rank);
    yp = vect(n);
    z = vect(max_rank);
    ryp = vect(max_rank);
    gptol = tolchol;
    cpar = vect(2);
    lpenMat = mymatrix(p, n);
    lpen = vect(n);
    activestatus = ivect(nsamp*p);
    
    v = vect(max_rank); u = vect(max_rank); mu_hat = vect(n + nnew); znorm = vect(max_rank); Hz = vect(max_rank), Rz = vect(n + nnew);
    GF = mymatrix(max_rank, max_rank); H = mymatrix(max_rank, max_rank); HC = mymatrix(max_rank, max_rank);
    
    int pos = 0, pos2 = 0;
    for(j = 0; j < p; j++){
        for(i = 0; i < n; i++) x[i][j] = xvar[pos + i]; pos += n;
        for(i = 0; i < nnew; i++) x[n + i][j] = xnew[pos2 + i]; pos2 += nnew;
    }
    
    pos = 0;
    for(j = 0; j < p; j++){
        for(i = 0; i < n; i++) lpenMat[j][i] = lpenalty[pos + i]; pos += n;
    }
    
    int *ix = ivect(p*ncomp), shift1 = 0, shift2 = 0, shift3 = 0;
    char b[10000];
    GetRNGstate();
    for(i = 0; i < nsamp; i++){
        for(l = 0; l < ncomp; l++){
            for(j = 0; j < p; j++) ix[l*p + j] = 0;
            strcpy(b, ix_store[shift1]);
            string_locator(b, ix + l*p, p);
            shift1++;
        }
        sigstore[shift3] = ragpFn(ncomp, 2, ix, y, par_store + shift2, fmeanstore + i*N, fstore + i*N); shift3++;
        shift2 += 2*ncomp;
    }
    PutRNGstate();
}
