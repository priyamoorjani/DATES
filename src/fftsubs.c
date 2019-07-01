#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <math.h>
#include <fftw3.h>

#include <nicklib.h>

/** 
 uses fftw3 library 
 3 main routines 
 fftcv Compute non-negative lag correlations 
 fftcv2 Compute correlations  + and - ;  Lag 0 is returned twice
 fftauto Compute non-negative lag autocorrelations  (similar to fftcv but slightly more efficient) 
 All routines use power of 2 fft and a slight efficiency gain might be possible here
*/

void fftauto(double *cout, double *a1, int n, int maxlag) ;
void fftcv(double *cout, double *a1, double *b1, int n, int maxlag) ;
void fftcv2(double *cout, double *cout2, double *a1, double *b1, int n, int maxlag) ;

double cnorm2(fftw_complex a)
{
 double y  ; 

 y = a[0]*a[0]  + a[1]*a[1] ;
 return y ; 

}

void cmul(fftw_complex c, fftw_complex a, fftw_complex b)
// complex multiply
{

  c[0] = a[0]*b[0] - a[1]*b[1] ;
  c[1] = a[0]*b[1] + a[1]*b[0] ;

}
void copyfft(fftw_complex *z1, fftw_complex *z2, int n) 
{
   int k ; 

  for (k=0; k<n; ++k) { 
     z2[k][0] = z1[k][0] ;
     z2[k][1] = z1[k][1] ;
  }

}

void slocv(double *cout, double *a1, double *b1, int m, int maxlag) 
{

  int i,j,k ;  
 
  vzero(cout, maxlag+1) ;  
  for (i=0; i<m; i++) { 
   for (k=0; k<=maxlag; ++k) { 
    j = i+k ; 
    if (j>=m) break ; 
    cout[k] += a1[i]*b1[j] ;
   }
  }
}

void fftcv(double *cout, double *a1, double *b1, int m, int maxlag) 
{

  fftw_complex *z1, *z2, *zz1, *zz2, *zin, *zout ;  
  fftw_plan p1, p2 ;
  int k, n, t ;
  double y1, y2, yn ;

  if (a1==b1) { 
    fftauto(cout, a1, m, maxlag) ;
    return ; 
  }

  y1 = log2(m) ;  
  y1 = ceil(y1)+1 ;  
  y2 = pow(2, y1) ; 
  n  = nnint(y2) ;
  yn = (double) n ;

  ZALLOC(z1, n, fftw_complex) ;
  ZALLOC(z2, n, fftw_complex) ;
  ZALLOC(zz1,n,  fftw_complex) ;
  ZALLOC(zz2,n,  fftw_complex) ;
  ZALLOC(zin,n,  fftw_complex) ;
  ZALLOC(zout,n,  fftw_complex) ;
 
  for (k=0; k<m; ++k) { 
   z1[k][0] = a1[k] ; 
   z2[n-k-1][0] = b1[k] ;
  }
 
  p1  = fftw_plan_dft_1d(n, zin, zout, FFTW_FORWARD, FFTW_ESTIMATE);
  p2  = fftw_plan_dft_1d(n, zin, zout, FFTW_BACKWARD, FFTW_ESTIMATE);

  if (p1==NULL) fatalx("fft trouble 1\n") ;
  if (p2==NULL) fatalx("fft trouble 2\n") ;

  copyfft(z1, zin, n) ;
  fftw_execute(p1) ;         
  copyfft(zout, zz1, n) ; 

  copyfft(z2, zin, n) ;
  fftw_execute(p1) ;         
  copyfft(zout, zz2, n) ; 

  for (k=0; k<n; ++k) { 
    cmul(z2[k], zz1[k], zz2[k]) ;
  }

  copyfft(z2, zin, n) ;
  fftw_execute(p2) ;

  vzero(cout, maxlag+1) ;
  for (k=0; k<=maxlag; ++k) { 
   t = n - k - 1 ; 
   if (t>=0)  cout[k] = zout[t][0] ;
  }

  vst(cout, cout, 1.0/yn,  maxlag+1) ;

  free(z1) ;
  free(z2) ;
  free(zz1) ; 
  free(zz2) ;
  free(zin) ; 
  free(zout) ;

  fftw_destroy_plan(p1) ;
  fftw_destroy_plan(p2) ;

}
void fftcv2(double *cout, double *cout2, double *a1, double *b1, int m, int maxlag) 
{

  fftw_complex *z1, *z2, *zz1, *zz2, *zin, *zout ;  
  fftw_plan p1, p2 ;
  int k, j, n, t ;
  double y1, y2, yn ;

  y1 = log2(m) ;  
  y1 = ceil(y1)+1 ;  
  y2 = pow(2, y1) ; 
  n  = nnint(y2) ;
  yn = (double) n ;

  ZALLOC(z1, n, fftw_complex) ;
  ZALLOC(z2, n, fftw_complex) ;
  ZALLOC(zz1,n,  fftw_complex) ;
  ZALLOC(zz2,n,  fftw_complex) ;
  ZALLOC(zin,n,  fftw_complex) ;
  ZALLOC(zout,n,  fftw_complex) ;
 
  for (k=0; k<m; ++k) { 
   z1[k][0] = a1[k] ; 
   z2[n-k-1][0] = b1[k] ;
  }
 
  p1  = fftw_plan_dft_1d(n, zin, zout, FFTW_FORWARD, FFTW_ESTIMATE);
  p2  = fftw_plan_dft_1d(n, zin, zout, FFTW_BACKWARD, FFTW_ESTIMATE);

  copyfft(z1, zin, n) ;
  fftw_execute(p1) ;         
  copyfft(zout, zz1, n) ; 

  copyfft(z2, zin, n) ;
  fftw_execute(p1) ;         
  copyfft(zout, zz2, n) ; 

  for (k=0; k<n; ++k) { 
    cmul(z2[k], zz1[k], zz2[k]) ;
  }


  copyfft(z2, zin, n) ;
  fftw_execute(p2) ;


  vzero(cout, maxlag+1) ; 
  vzero(cout2, maxlag+1) ; 

  for (k=0; k<=maxlag; ++k) { 
   t = n - k - 1 ; 
   if (t>=0) cout[k] = zout[t][0] ;
   j = k-1 ; 
   if (j<0) j = n-1 ; 
   if (j<n) cout2[k] = zout[j][0] ;
  }
  vst(cout, cout, 1.0/yn,  maxlag+1) ;
  vst(cout2, cout2, 1.0/yn,  maxlag+1) ;

  free(z1) ;
  free(z2) ;
  free(zz1) ; 
  free(zz2) ;
  free(zin) ; 
  free(zout) ;

  fftw_destroy_plan(p1) ;
  fftw_destroy_plan(p2) ;

}
void fftauto(double *cout, double *a1,  int m, int maxlag) 
{

  fftw_complex *z1, *z2, *zz1, *zz2, *zin, *zout ;  
  fftw_plan p1, p2 ;
  int k, n, nn, outmax ;
  double y1, y2, yn ;

  y1 = log2(m) ;  
  y1 = ceil(y1)+1 ;  
  y2 = pow(2, y1) ; 
  n  = nnint(y2) ;
  yn = (double) n ;
//  printf("zz %d\n", n) ;
  nn = MAX(n, maxlag+1) ;

  ZALLOC(z1, nn, fftw_complex) ;
  ZALLOC(z2, nn, fftw_complex) ;
  ZALLOC(zz1,nn,  fftw_complex) ;
  ZALLOC(zz2,nn,  fftw_complex) ;
  ZALLOC(zin,nn,  fftw_complex) ;
  ZALLOC(zout,nn,  fftw_complex) ;
 
  for (k=0; k<m; ++k) { 
   z1[k][0] = a1[k] ; 
  }
 
  p1  = fftw_plan_dft_1d(n, zin, zout, FFTW_FORWARD, FFTW_ESTIMATE);
  p2  = fftw_plan_dft_1d(n, zin, zout, FFTW_BACKWARD, FFTW_ESTIMATE);

  copyfft(z1, zin, n) ;
  fftw_execute(p1) ;         
  copyfft(zout, zz1, n) ; 

  for (k=0; k<n; ++k) { 
    z2[k][0] = cnorm2(zz1[k]) ; 
    z2[k][1] = 0 ; 
  }

  copyfft(z2, zin, n) ;
  fftw_execute(p2) ;

  vzero(cout, maxlag+1) ; 

  outmax = MIN(n-1, maxlag) ;
  for (k=0; k<=outmax; ++k) { 
   cout[k] = zout[k][0] ;
  }

  vst(cout, cout, 1.0/yn,  outmax) ;

  free(z1) ;
  free(z2) ;
  free(zz1) ; 
  free(zz2) ;
  free(zin) ; 
  free(zout) ;

  fftw_destroy_plan(p1) ;
  fftw_destroy_plan(p2) ;

}
