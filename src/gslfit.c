#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <nicklib.h>
#include <getpars.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_errno.h>
#include "qpsubs.h"

#define WVERSION   "100"
#define MAXSTR  512

extern int gsldetails;
extern double gslprecision;

double scorexp (double *www) ;
static int xnexp;
static double *xvexp;
static double *www;

gsl_multimin_fminimizer *s = NULL;
static gsl_vector *ss, *x;
gsl_multimin_function minex_func;

double fjunk2 (const gsl_vector *v, void *params);
static double fff (double *a, int n);

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


int
gslsetup (int nexp, double *vexp)
{
  int k;
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;

  if (nexp == 0)
    return 0;
  xnexp = nexp;
  ZALLOC (xvexp, nexp, double);
  copyarr (vexp, xvexp, nexp);
  ZALLOC (www,  nexp, double);


  minex_func.n = nexp;
  minex_func.f = fjunk2;
  minex_func.params = NULL;

  x = gsl_vector_alloc (nexp);
  for (k = 0; k < nexp; ++k) {
    gsl_vector_set (x, k, s2p(vexp[k]));
  }

  /* Set initial step sizes to 0.01 */
  ss = gsl_vector_alloc (nexp);
  gsl_vector_set_all (ss, 0.01);

  /* Initialize method and iterate */

  s = gsl_multimin_fminimizer_alloc (T, nexp);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  gsl_set_error_handler_off ();

//  gsldetails = YES ;  
  printf ("gslsetup called\n");
  
  return 1;
}


double
gslopt (double *wpars)
{
  size_t iter = 0, k;
  int status;
  double size;
  double q = 999999;
  double qbest = 1.0e40;

  /* Starting point */


  for (k = 0; k < xnexp; ++k) {
    gsl_vector_set (x, k, s2p(wpars[k]));
  }
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  gsl_set_error_handler_off ();

  do {
    iter++;
    status = gsl_multimin_fminimizer_iterate (s);

    if (status)
      break;

    size = gsl_multimin_fminimizer_size (s);
    status = gsl_multimin_test_size (size, gslprecision);
    q = s->fval;
    if (gsldetails)
      printf ("gslopt: %3d %12.6f\n", (int) iter, q);
    if (q < qbest) {
      if (gsldetails)
	printf ("+++ new best\n");
      qbest = q;
      for (k = 0; k < xnexp; ++k) {
	wpars[k] = gsl_vector_get (s->x, k);
      }
    }
  }
  while (status == GSL_CONTINUE && iter < 10000);

  printf ("gslans: %4d %12.6f\n", (int) iter, q);
  fff (wpars, xnexp);
  printmat (wpars, 1, xnexp);

  fflush (stdout);
  return q;
}

double
fff (double *a, int n)
// make a cononical; return penalty
{
  double penalty = 0;
  int k;
  double xx, y, yy;

  for (k = 0; k < n; ++k) {
    yy = xx = a[k];
    if ((xx >= 0) && (xx <= 1))
      continue;
    if (xx < 0) {
      yy = -yy;
    }
    if (yy > 1) {
      yy = 2.0 * modf (0.5 * yy, &y);   // periodicity is 2
      if (yy > 1)
        yy = 2 - yy;
    }
    y = yy - xx;
    penalty += y * y;
    a[k] = xx = yy;
  }
  return penalty;
}


double
fjunk2 (const gsl_vector * v, void *params)
{
  double xx;
  double q;
  int k, t = 0;
  double *tt;
  double penalty = 0;

  ZALLOC (tt, xnexp, double);
  for (k = 0; k < xnexp; ++k) {
    tt[k] = gsl_vector_get (v, k);
  }
  penalty = fff (tt, xnexp);

  t = 0;
  for (k = 0; k < xnexp; ++k) {
    xx = tt[k];
    www[k] = xx ;     
  }
  q = scorexp(www) ;
  q += 10 * penalty;
  if (gsldetails) {
    printf ("zzfjunk2 %12.6f\n", q);
    printmat (tt, 1, xnexp);
  }
  free (tt);
  return q;
}
