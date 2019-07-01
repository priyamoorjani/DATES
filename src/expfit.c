#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <nicklib.h>
#include <getpars.h>

#include "badpairs.h"
#include "admutils.h"
#include "mcio.h"  
#include "mcmcpars.h"  
#include "egsubs.h"  
#include "qpsubs.h"  

#define WVERSION   "200" 

/** 
 fit exponential
 gsl version
*/

double gslprecision = .000001 ;
int gsldetails = NO ;

#define MAXFL  50   
#define MAXSTR  512
#define MAXPOPS  100
#define NEXP  2  

int debug = NO ; 

extern int packmode ;

int popcheck = YES ;
int fancynorm = YES ; 
int plotmode = NO ;
int outnum = -1 ;
int affmode = NO ;
int helpmode = NO ;
int pop2mode = NO  ;
int popallmode = NO  ;
int numexp = NEXP  ;
char *iname = NULL ;
char *oname = NULL ;

int datacol = 1 ;
double loval = -1.0e20 ;
double hival = 1.0e20 ;
double addx  = 0 ;
double step = -1 ;

int popsizelimit = -1 ;
int xscratch[50000] ; 
double yscratch[100000] ;

char *trashdir = "/var/tmp" ;
int verbose = NO ;
int qtmode = NO ;
Indiv **indivmarkers;
SNP **snpmarkers ;
int numsnps, numindivs ; 

char  *genotypename = NULL ;
char  *snpname = NULL ;
char  *genooutfilename = NULL ;
char  *indoutfilename = NULL ;
char  *indivname = NULL ;
char *badsnpname = NULL ;
char *goodsnpname = NULL ;
char *badpairsname = NULL ;
char *markername = NULL ;
char *poplistname = NULL ;

char *outputname = NULL ;
FILE *ofile ;

int xchrom = -1 ;
int *fchrom, *lchrom ;

int details = NO ;
double fakespacing = 0.0 ;
double dthresh =  0.02 ;
int    heminum  = 50 ;
int    homonum  = 50 ;
double blgsize  = .001 ;  // morgans
double blglim   = .10  ;  // morgans
long seed = 0 ;

double *xfit, *xexp, *xco, *xval, *xbase ;
double *zfit, *zexp, *zans ;
int nexp, lenxfit ;

char  unknowngender = 'U' ;

void readcommands(int argc, char **argv) ;
int dohemiz(Indiv *ind1, Indiv *ind2, SNP **snpmarkers, int numsnps, SNP ***runlist, int *lrunlist, double dthresh)  ;
void pubrun(SNP ***runlist, SNP *slo, SNP *shi, int *pnruns, double dthresh) ;
void printruns(SNP ***runlist, int nruns, Indiv *ind1, Indiv *ind2) ;
double fitexp(double *xval, double *xfit, double *xexp, double *xco, int m, int n, int mode, int inititer) ;

void dofit(double **xx, int n)  ;
void printhelp() ;
double scorit(double *x, double *vv, int m, int n, int mode, double *zans, double *zexp, double *zfit) ;
double scorexp(double *vexp) ;
double p2s(double x) ;
double s2p(double x) ;

int main(int argc, char **argv)
{

  int i, j, k, g, t, ipop, jpop, kpop ; 
  int k1, k2 ;
  int len ;

  int ch1, ch2 ;

  int numvind, nignore, numrisks = 1 ;
  int numcols ;
  int maxbin ;
  double y, gdis, ysd ;
  double **xx ;  
  char *pop1, *pop2 ;
  char popstring[MAXSTR] ;   
  int startiter ; 


  ofile = stdout; 
  readcommands(argc, argv) ;
  printf("%s version: %s\n", argv[0], WVERSION) ;
  if (helpmode) {
    printhelp() ;
    return 1 ; 
  } 

  if (iname==NULL) {
     printhelp() ;
     return 0 ;
  }
  if (oname != NULL) openit(oname, &ofile, "w") ;

  nexp = numexp ;
  maxbin = numlines(iname) ; 

  numcols = datacol+1 ; 
  xx = initarray_2Ddouble(numcols, maxbin,  0.0) ;

  maxbin = getxx(xx, maxbin, numcols, iname) ;

  xbase = xx[0] ;
  xval = xx[datacol] ;
  lenxfit = maxbin ; 
  step = (xbase[1]-xbase[0])/100.0  ;  // step in Morgans
  printf("step (Morgans) :: %12.6f\n", step) ;
  if (step < 0) fatalx("step negative!\n") ; 
  for (;;) { 
   y = xbase[0] ; 
   if (y>=loval) break ;
   ++xbase ; 
   ++xval ; 
   --lenxfit ;
   if (lenxfit <= 0) fatalx("no data\n") ;
  }

  len = 0 ;
  for (k=0; k<lenxfit; ++k) {
   y = xbase[k] ; 
   if (y>hival) break ;
   ++len ;
  }
  lenxfit = len ;

  vsp(xval, xval, addx, lenxfit) ;
  ZALLOC(xfit, lenxfit, double) ;
  ZALLOC(xexp, nexp+1, double) ;
  ZALLOC(xco, nexp+1, double) ;

  ZALLOC(zexp, nexp+1, double) ;  
  ZALLOC(zans, nexp+1, double) ;  
  ZALLOC(zfit, lenxfit+1000, double) ;  

  k = numexp ;
   if (affmode) printf("fitting %d exponentials + affine\n", k) ;
   if (!affmode) printf("fitting %d exponentials\n", k) ;
    y = 100*pow(2, numexp) ; 
    startiter = nnint(y) ;
    ysd = fitexp(xval, xfit, xexp, xco, k, lenxfit, affmode, startiter) ;
    fflush(stdout) ; 
    y = scorexp(xexp) ; 

//  printf("zz1 %12.6f %12.6f\n", ysd, sqrt(y)) ;
/**
     printf("exp, co nexp: %d\n", nexp ) ;
     printmatl(xexp, 1, nexp) ;
     printmatl(xco, 1, nexp+1) ;
*/
    
    for (i=0; i<nexp; i++) { 
     xexp[i] = hlife(xexp[i]) ; 
    }
    printf("error sd: %12.6f\n", ysd) ;
    printf("halflife: ") ;
    printmat(xexp, 1, nexp) ;
    if (step>0) {  
     for (i=0; i<nexp; i++) { 
      y = xexp[i]*step;   
      y += 1.0e-20 ;
      xexp[i] = log(2.0) / y ;
     }
     printf("mean (generations): ") ;
     printmat(xexp, 1, nexp) ;
     for (i=0; i<nexp; ++i) { 
      xco[i] *= exp(xexp[i]*(xbase[0]*.01-step)) ;
     }
     printmatl(xco, 1, nexp+1) ;
    }
    fprintf(ofile, "##fit: %s\n", iname) ;
    for (i=0; i<lenxfit; ++i)  { 
     fprintf(ofile, "%12.6f ", xbase[i]) ;
     fprintf(ofile, "%12.6f ", xval[i]) ;
     fprintf(ofile, "%12.6f ", xfit[i]) ;
     fprintf(ofile, "%12.6f ", xval[i] - xfit[i]) ;
     fprintf(ofile, "\n") ;
    }
  printf("##end of run\n") ;
  return 0 ;
}
void printhelp() 
{ 
 printf("expfit: \n") ;
 printf (" -i iname   ## input\n") ;
 printf (" -o oname   ## output\n") ;
 printf (" -n numexp  ## number of exponentials\n") ;
 printf (" -c col     ## data column  (0 is xval)\n") ;
 printf (" -l loval   ## lowest x value\n") ;
 printf (" -h hival   ## highest x value\n") ;
 printf (" -s stepsize   ## step size (Morgans)\n") ;
 printf (" -x val     ## value to add to x (deprecated)\n") ;
 printf (" -r ran     ## seed for random generator\n")  ;
 printf (" -a         ## affine mode (add constant)\n") ;
 printf (" -V         ## verbose mode \n") ;
 printf (" -m         ## print help menu and quit\n") ;

}


void readcommands(int argc, char **argv) 

{
  int i,haploid=0;
  char *parname = NULL ;
  int n ;

  while ((i = getopt (argc, argv, "c:i:n:r:s:o:l:h:x:Vam")) != -1) {

    switch (i)
      {

      case 'i':
	iname = strdup(optarg) ;
	break;

      case 'n':
	numexp = atoi(optarg) ;
	break;

      case 'r':
	seed = atoi(optarg) ;
	break;

      case 'c':
	datacol = atoi(optarg) ;
	break;

      case 's':
	step = atof(optarg) ; // morgans
	break;

      case 'a':
	affmode = YES ;
	break; 

      case 'm':
	helpmode = YES ;
	break; 

      case 'o':
	oname = strdup(optarg) ;
	break; 

      case 'l':
	loval = atof(optarg) ;
	break; 

      case 'h':
	hival = atof(optarg) ;
	break; 

      case 'x':
	addx = atof(optarg) ;
	break; 

      case 'V':
	verbose = YES ;
	break; 

      default: 
	printf ("Usage: bad params.... \n") ;
	fatalx("bad params\n") ;
      }
  }
         
}

double scorexp(double *vexp) 
{

  double y ; 
  double *params ; 
  int i ;  
  
  y = scorit(xval, vexp, numexp, lenxfit, affmode, zans, zexp, zfit) ; 
  return y ;

  ZALLOC(params, nexp, double) ;
  for (i=0; i<numexp; i++) { 
   params[i] = s2p(vexp[i]) ;
  }
  y = scorit(xval, params, numexp, lenxfit, affmode, zans, zexp, zfit) ; 
  free(params) ;
  return y ;

}

double scorit(double *x, double *vv, int m, int n, int mode, double *zans, double *zexp, double *zfit) 
{
   int i, mm ; 
   double *eq, *rhs ;
   double z, y1, y2 ; 
   double *tt ;

   mm = m ; 
   if (mode) ++mm ;
   ZALLOC(eq, mm*n, double) ;
   ZALLOC(rhs, n, double) ;
   ZALLOC(tt, mm, double) ;

   if (debug) { 
    printf("zzdebug: ");  printmat(vv, 1, m) ; 
   }

   vclear(tt, 1.0, mm) ;
   for (i=0; i<n; i++) { 
    vvt(tt, tt, vv, m) ;
    copyarr(tt, eq+i*mm, mm) ;
    rhs[i] = x[i] ;
   }


   regressit(tt, eq, rhs, n, mm) ;
   y2 = 0.0  ;                            
   for (i=0; i<n; i++)  {
    zfit[i] = vdot(tt, eq+i*mm, mm) ;  
    z = rhs[i] - zfit[i] ;
    y2 += z*z ;
   }
   y2 /= (double) n ;
 
   free(eq) ;
   free(rhs) ;

   copyarr(tt, zans, mm) ;
   copyarr(vv, zexp, m) ;
   free(tt) ;

   if (debug) { 
    printf("zzd2 %15.9f ", y2) ;
    printmat(vv, 1, m) ;
   }

   return y2 ;

}
