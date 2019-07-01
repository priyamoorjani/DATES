#include <stdio.h>
#include <limits.h>
#include <math.h>  
#include <nicklib.h>

extern int verbose ;

extern double *zans, *zfit, *zexp ;
static double  *xxval = NULL ; 
static int  mval, nval, llag, affmode = YES ;

double scorit(double *x, double *vv, int m, int n, int mode, double *zans, double *zexp, double *zfit) ;
double p2s (double p) ;
double s2p (double s) ;

void loadval(double *x, int m, int n, int nlag, int mode) ;

int gslsetup (int nexp, double *vexp) ;
double gslopt (double *wpars) ;

extern long seed  ; 
extern int debug ; 

double fitexp(double *xval, double *xfit, double *xexp, double *xco, int m, int n, int mode, int inititer) 
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

  ZALLOC(try, m, double) ;
  ZALLOC(tbest, m, double) ;
  ZALLOC(ww, m, double) ;
  ZALLOC(wbase, m, double) ;

  vbest = 1.0e6 ; 
  loadval(xval, m, n, 0, mode) ;
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
   val = scorit(xval, try, m, n, affmode, zans, zexp, zfit) ;
    if (val < vbest)  { 
      vbest = val ;
      copyarr(try, tbest, m) ;
//    printf("%12.6f \n", vbest*1.0e6) ;
//    printmat(try, 1, m) ;
    }
    x = iter % inniter ;  
   if (x == 0) mul *= 0.5 ;
   
   }
  printf("after initialization: %12.6f ", vbest);  printmat(tbest, 1, m) ;
  verbose = NO ;
  
  gslsetup (m, tbest) ; 
  gslopt(tbest) ;

 // debug = YES ;  
  val = scorit(xval, tbest, m, n, affmode, zans, zexp, zfit) ;
 if (debug) {
  printf("zzvvv %12.6f \n", sqrt(val)) ;
  printmatl(tbest, 1, m) ;
 }
  mm = m ; 
  if (affmode) ++mm ;
  
  copyarr(zexp, xexp, m) ;
  copyarr(zans, xco, mm) ;
  copyarr(zfit, xfit, n) ;
  free(try) ;
  free(tbest) ;
  free(ww) ;
  free(wbase) ;
  val += 1.0e-12 ;
  return sqrt(val) ;
}

void loadval(double *x, int m, int n, int nlag, int mode) 
{
  if (xxval != NULL) free(xxval)  ;
  ZALLOC(xxval, n, double) ;
  mval = m ; 
  nval = n ;
  llag = nlag ;
  affmode = mode ;

  copyarr(x, xxval, n) ;

}



