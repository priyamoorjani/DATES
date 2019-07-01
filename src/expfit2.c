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

#define WVERSION   "110" 

// badpairsname added 
/** 
 fit exponential

 New I/O (mcio.c) added
*/


#define MAXFL  50   
#define MAXSTR  512
#define MAXPOPS  100
#define NEXP  2  


extern int packmode ;

int popcheck = YES ;
int fancynorm = YES ; 
int plotmode = NO ;
int outnum = -1 ;
int affmode = NO ;
int pop2mode = NO  ;
int popallmode = NO  ;
int numexp = NEXP  ;
char *iname = NULL ;
char *oname = NULL ;

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

char  unknowngender = 'U' ;

void readcommands(int argc, char **argv) ;
int dohemiz(Indiv *ind1, Indiv *ind2, SNP **snpmarkers, int numsnps, SNP ***runlist, int *lrunlist, double dthresh)  ;
void pubrun(SNP ***runlist, SNP *slo, SNP *shi, int *pnruns, double dthresh) ;
void printruns(SNP ***runlist, int nruns, Indiv *ind1, Indiv *ind2) ;
void fitexp2(double *xval, double *xfit, double *xexp, double *xco, 
   int m, int n, int mode, int inititer) ;

void dofit(double **xx, int n)  ;

int main(int argc, char **argv)
{

  int i, j, k, g, t, ipop, jpop, kpop ; 
  int k1, k2, len ;

  int ch1, ch2 ;

  int numvind, nignore, numrisks = 1 ;
  int numcols ;
  int maxbin ;
  double y, gdis ;
  double **xx ;  
  char *pop1, *pop2 ;
  char popstring[MAXSTR] ;   

  double *xfit, *xexp, *xco, *xval, *xbase ;
  double yaff ;
  int nexp, lenxfit ;

  ofile = stdout; 
  readcommands(argc, argv) ;
  if (oname != NULL) openit(oname, &ofile, "w") ;

  nexp = numexp ;
  if (iname==NULL) fatalx("i param compulsory") ;
  maxbin = numlines(iname) ; 

  numcols = 2 ; 
  xx = initarray_2Ddouble(numcols, maxbin,  0.0) ;

  maxbin = getxx(xx, maxbin, numcols, iname) ;

  xbase = xx[0] ;
  xval = xx[1] ;
  lenxfit = maxbin ; 
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
  ZALLOC(xco, 2*nexp+1, double) ;
  k = numexp ;
   if (affmode) printf("fitting %d exponentials (d 2d) + affine\n", k) ;
   if (!affmode) printf("fitting %d exponentials (d 2d) \n", k) ;
   fitexp2(xval, xfit, xexp, xco, k, lenxfit, affmode, 5000) ;
    if (affmode) {  
     yaff = xco[2*nexp] ; 
     printf("affine: %12.6f\n", yaff) ;
    }
    printf("exp, co nexp: %d\n", nexp ) ;
    printmatl(xexp, 1, nexp) ;
    printmatl(xco, 1, nexp) ;
    printmatl(xco+nexp, 1, nexp) ;
    for (i=0; i<nexp; i++) { 
     xexp[i] = hlife(xexp[i]) ; 
    }
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

void readcommands(int argc, char **argv) 

{
  int i,haploid=0;
  char *parname = NULL ;
  int n ;

  while ((i = getopt (argc, argv, "i:n:r:s:o:l:h:x:a")) != -1) {

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

      case 's':
	step = atof(optarg) ; // morgans
	break;

      case 'a':
	affmode = YES ;
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

      case '?':
	printf ("Usage: bad params.... \n") ;
	fatalx("bad params\n") ;
      }
  }

         
}

void dophyscheck(SNP **snpm, int numsnps) 
{
// catch places where physpos genpos are in opposite order
  SNP *cupt, *cuptold ;
  int i ;

  for (i=0; i<numsnps; i++) {   
   cupt = snpm[i] ;
   if (i==0) cuptold = cupt ;
   if (cupt -> isfake) continue ;
   if (cupt -> ignore) continue ;
   if (cupt -> chrom == cuptold -> chrom)  {
    if (cupt -> physpos < cuptold -> physpos) {  
     printf("physcheck %20s %15s %12.3f %12.3f %13.0f %13.0f\n", 
     cuptold->ID, cupt -> ID, 
     cuptold -> genpos, cupt -> genpos, 
     cuptold -> physpos, cupt -> physpos);
    }
   }
   cuptold = cupt ;
  }
}
void cleartagnumber(SNP **snpm, int numsnps) 
{
  int i ; 
  for (i=0; i<numsnps; i++) {   
   snpmarkers[i] -> tagnumber = 0 ;
  }
}

void setfc(SNP **snpm, int numsnps) 
// also allocates fchrom , lchrom
{
  int i, j, chrom ;
  SNP *cupt ;
  double dis, totdis, pos1, pos2 ;

  ZALLOC(fchrom, 25, int) ;
  ZALLOC(lchrom, 25, int) ;

  ivclear(fchrom, 999999, 25) ;
  ivclear(lchrom, -999999, 25) ;
  
/* initialize real marker array */
  for (i=0; i<numsnps; i++) {
    cupt = snpm[i] ;
//  if (cupt -> isfake) continue ;
    if (cupt -> ignore) continue ;
    chrom = cupt -> chrom ;
    fchrom[chrom] = MIN(fchrom[chrom],i) ;
    lchrom[chrom] = MAX(lchrom[chrom],i) ;
  }

  totdis = 0.0 ;
   if(details)
     {
       printf("\n###GENETIC DISTANCE FOR ALL CHROMOSOMES\n");
       printf("##Chr_Num: chromosome num, First_SNP and Last_SNP: First and last markers, Gen_dist: Genetic distance\n");
       printf("%5s %9s %9s %9s\n","Chr_Num","First_SNP","Last_SNP","Gen_dist");
     }

  for (j=1; j<=23; j++) {  
   if (lchrom[j]<0) continue ;
   cupt = snpm[fchrom[j]]; 
   pos1 = cupt -> genpos ;
   cupt = snpm[lchrom[j]]; 
   pos2 = cupt -> genpos ;
   dis = pos2 - pos1 ;
   totdis += dis ;
   printf("chrom: %5d  first: %5d  last:  %3d", 
     j, fchrom[j], lchrom[j]) ;
   printf(" distance: %9.3f\n", dis) ;
  }
  printf("total distance: %9.3f\n", totdis) ;
}

