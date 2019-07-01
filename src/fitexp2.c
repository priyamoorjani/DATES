#include <nag.h>
#include <nag_stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>  
#include <nicklib.h>
#include <nage04.h>
#include "regsubs.h" 

extern int verbose ;

double *zans, *zfit, *zexp ;
static double  *xxval = NULL ; 
static int  mval, nval, affmode = YES ;

double scorit(double *x, double *vv, int m, int n, int mode) ;
double p2s (double p) ;
double s2p (double s) ;
void loadval(double *x, int m, int n, int mode) ;
void nagopt(double *params, int n) ;
extern long seed  ; 

double fitexp2(double *xval, double *xfit, double *xexp, double *xco, 
  int m, int n, int mode, int inititer) 
{

  int i, k, numcols, x, iter, mm, outiter = 10, inniter   ;
  double y, ytry, val, vbest, ybest, yinc ;
  char *fname = "qq" ;
  double **xx, **fit, *try, *tbest, *ww, *wbase ;
  double lo = 0.0, mul = 1.0  ;

  if (seed == 0) {
   seed = seednum() ;
   printf("seed: %ld ", seed) ; printnl() ;
   SRAND(seed)  ;
  }

  mm = 2*m ; 
  if (affmode) ++mm ;

  ZALLOC(try, m, double) ;
  ZALLOC(tbest, m, double) ;
  ZALLOC(zexp, m, double) ;
  ZALLOC(zans, 2*m+1, double) ;
  ZALLOC(zfit, n, double) ;
  ZALLOC(ww, m, double) ;
  ZALLOC(wbase, m, double) ;

  vbest = 1.0e6 ; 
  loadval(xval, m, n, mode) ;
  inniter = inititer / outiter ;
  inniter += outiter ;
  for (iter = 1; iter <= inititer; ++iter) {
    for (k=0; k<m; ++k) {  
     try[k] = lo + (1.0-lo) *DRAND() ;
    }

   vst(ww, tbest, 1.0-mul, m) ;  
   vst(try, try, mul, m) ;
   vvp(try, try, ww, m) ;

   sortit(try, NULL, m) ;
   val = scorit(xval, try, m, n, affmode) ;
    if (val < vbest)  { 
      vbest = val ;
      copyarr(try, tbest, m) ;
      printf("%12.6f \n", vbest*1.0e6) ;
      printmat(try, 1, m) ;
    }
    x = iter % inniter ;  
   if (x == 0) mul *= 0.5 ;
   
   }
  printmat(tbest, 1, m) ;
  verbose = NO ;
  nagopt(tbest, m) ;
  verbose = NO ;
  val = scorit(xval, tbest, m, n, affmode) ;
  printf("%12.6f \n", val*1.0e6) ;
  printmatl(tbest, 1, m) ;
  
  copyarr(zexp, xexp, m) ;
  copyarr(zans, xco, mm) ;
  copyarr(zfit, xfit, n) ;
  free(try) ;
  free(tbest) ;
  free(zexp) ;
  free(zans) ;
  free(zfit) ;
  free(ww) ;
  free(wbase) ;
  val += 1.0e-12 ; 
  return sqrt(val) ;
}

double scorit(double *x, double *vv, int m, int n, int mode) 
{
   int i, mm ; 
   double *eq, *rhs, *v2 ;
   double z, y1, y2 ; 
   double *tt ;

   mm = 2*m ; 
   if (mode) ++mm ;
   ZALLOC(eq, mm*n, double) ;
   ZALLOC(rhs, n, double) ;
   ZALLOC(tt, mm, double) ;
   ZALLOC(v2, mm, double) ;
   copyarr(vv, v2, m) ;
   vvt(v2+m, vv, vv, m) ;

   vclear(tt, 1.0, mm) ;
   for (i=0; i<n; i++) { 
    vvt(tt, tt, v2, 2*m) ;
    copyarr(tt, eq+i*mm, mm) ;
    rhs[i] = x[i] ;
   }

   regressit(tt, eq, rhs, n, mm) ;
   copyarr(tt, zans, mm) ;
   copyarr(vv, zexp, m) ;
   y2 = 0.0  ;                            
   for (i=0; i<n; i++)  {
    zfit[i] = vdot(tt, eq+i*mm, mm) ;  
    z = rhs[i] - zfit[i] ;
    y2 += z*z ;
   }
   y2 /= (double) n ;
 
   free(eq) ;
   free(rhs) ;
   free(tt) ;
   free(v2) ;

   if (verbose) { 
    printf("zz %15.9f ", y2) ;
    printmat(vv, 1, m) ;
   }
   return y2 ;

}
double p2s (double p) 
{
  if (p==0.0) return -1.0e20 ;
  if (p==1.0) return 1.0e20 ;
  return log(p/(1-p)) ;
}

double s2p (double s) 
{
  double r ;
  r = exp(s) ; 
  return  r/(1+r) ;
}


void loadval(double *x, int m, int n, int mode) 
{
  if (xxval != NULL) free(xxval)  ;
  ZALLOC(xxval, n, double) ;
  mval = m ; 
  nval = n ;
  affmode = mode ;

  copyarr(x, xxval, n) ;


}

void myfun(Integer n, double *xc, double *fc, Nag_Comm *comm) 
{  
   int nparams = n ;
   double *params ; 
   double y, yval, yfix ;
   int i  ;

   ZALLOC(params, mval, double) ;
   for (i=0; i<mval; i++) { 
    params[i] = s2p(xc[i]) ;
   }
   yval = 1.0e6*scorit(xxval, params, mval, nval, affmode) ;                                                     
   free(params) ;
   *fc =  yval ;

}

void nagopt(double *params, int n) 
{

  double val ;
  Nag_E04_Opt options;
  static NagError fail;
  double *ycoeffs ;
  static int ncall = 0 ; 
  int nparams = n, i ;
  double yval, oldyval, tmp ;


  ++ncall ;
  e04xxc(&options);

  options.print_level = Nag_NoPrint ;
  if (verbose) options.print_level = Nag_Soln ;
  options.optim_tol = 1.0e-10 ; 
  options.max_iter = n*2000 ;

  options.list = FALSE ; 

  fail.print = TRUE;

  ZALLOC(ycoeffs, nparams, double) ;
  for (i=0; i<mval; i++) { 
    ycoeffs[i] = p2s(params[i]) ;
  }

  oldyval = yval = scorit(xxval, params, mval, nval, affmode) ;                                           

/**
  printf("nagopt ncall: %d  value:  %9.3f\n", ncall, yval) ;

  for (i=0; i<n; i++) { 
   printf("%3d %12.6fn", i, params[i]) ;  
  }
*/

  e04ccc(n, myfun, ycoeffs, &yval, &options, NAGCOMM_NULL, &fail);
  for (i=0; i<mval; i++) { 
    params[i] = s2p(ycoeffs[i]) ;
  }
  yval = scorit(xxval, params, mval, nval, affmode) ;                                           
  printf("nagopt ncall: %d  oldvalue:  %15.9f newvalue: %15.9f\n", ncall, oldyval*1.0e6, yval*1.0e6) ;
  free(ycoeffs) ;
}




