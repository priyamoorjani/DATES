
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include  <math.h>

#include <nicklib.h>
#include <getpars.h>
#include <globals.h> 

#include "badpairs.h"
#include "admutils.h"
#include "mcio.h"  
#include "mcmcpars.h"  
#include "egsubs.h"  
#include "ldsubs.h"  

#define WVERSION   "753" 

// mincount stripped 
// basic idea: Model each admixed sample as a linear fit 
// of 2 parental pops.  Compute residual and rolloff (res*wt) 
// jackknife by chromosome
// weights added as option.  Latest I/O (mcio.c)
// timeoffset added 
// mesh + FFT convolutions.  slightly faster FFT
// zero popsize trapped.  runfit option
// admixlist introduced  
// snpname, seed passed into rolljack
// fft bug fixed (very low coverage) 
// PM: modified rolljack call to jackknife
// covariance output in col 1 of output file (dumpit) 
// trap zero snps in dowork 

#define MAXFL  50   
#define MAXSTR  512
#define MAXAPOPS  100

#define BINSIZE  .0005
#define MAXDIS   0.05  

extern int packmode ;
int zdipcorrmode = NO ; 

char *trashdir = "/var/tmp" ;
int qtmode = NO ;
int debug = NO ;
Indiv **indivmarkers;
SNP **snpmarkers ;
SNP **snpm ;
int numsnps, numindivs ; 
int checkmap = YES ;

char *parname = NULL ;

char  *genotypename = NULL ;
char  *snpname = NULL ;
char *snppath = NULL ; 
char  *genooutfilename = NULL ;
char  *indoutfilename = NULL ;
char  *indivname = NULL ;
char *badsnpname = NULL ;
char *poplistname = NULL ;
char *admixpop = NULL ;
char *admixlist = NULL ;
char **admixpoplist = NULL  ;
int nadmixpop = 0; 
char *dumpname = NULL ;
char *loadname = NULL ;
char *oname = NULL ;
char *weightname = NULL ; 
char *timeoffsetname = NULL ; 
int crossmode = NO ;
int runmode = 0 ;
int ldmode = NO ;
int ransample = -1 ;
int flatmode = NO  ;
double chithresh = -1.0; // default 6.0 
int jackknife = YES ;
int qbin = 0 ;

FILE *ofile ;
int  *fchrom, *lchrom ;
int  numeg = 0 ;
int nums = 2,   ns2, ns4 ;
int numbigiter = 2 ;
int xchrom = -1 ;
int xnochrom = -1 ;
int seed = 0 ;

int noxdata = YES ;
int samelambda = NO ;
int fixlambda = NO ;
int flatinit = NO ;  // different default
int popestimate = NO ; 
int nxlim = 100000 ;
int minparentcount = 10 ;
int norun = NO ;

double maxdis = MAXDIS   ;
double binsize = BINSIZE ; 

// if YES same across pops
int *firstsnp, *lastsnp ; 
int *firstqb, *lastqb ; 
int *qindex ;  

// expfit parameters 
int afffit = YES ; 
double lovalfit = .45 ;   
char *logfit = NULL ;
int runfit = YES ;


char  unknowngender = 'U' ;

int readcommands(int argc, char **argv) ;
void setfc(SNP **snpm, int numsnps) ;
void sethmm(double *theta, double lambda, int ind) ;
int addr4(int a, int b, int c, int d) ;
double doalpha(int ind) ;
double dobeta(int ind) ;

void cleargcounts()   ;
void reest()   ;
void calcgcounts(int ind)  ;
void printpars()   ;
void printgams()   ;
void calcgtl(int ind) ;
int  setqbins(SNP **snpm, int numsnps, double binsize, int qbin)  ; 
double doscore(double **thet, double *lam) ;
void dumppars(char *dumpname,  double **wmix, double **theta, int nums, int numeg) ;
void dorc(double *ans, double *res, double binsize, double maxdis) ;
void dorc2(double *ans, double *res, double binsize, double maxdis) ;
void dorcx(double *ans, double *rr, double *dd, double binsize, double maxdis) ;
void dorcxx(double *ans, double *rr, double *dd, double binsize, double maxdis) ;
void dumpit(char *oname, double *ww, CORR **corrbins, int len, double bsize, char *hdr) ;
double cntit1(double *xc, SNP *cupt, Indiv **indm, int numindivs, int t) ; 
double cntit2(double *xc, SNP *cupt, SNP *cupt2, Indiv **indm, int numindivs, int t) ; 
double estyran(SNP **snpmarkers, int numsnps, int nxlim) ;
void cx1(double *xc, double *y1, double *y2) ;
double cminor(double *xx) ;
int calcldscore(SNP *cupt, SNP *cupt2, double *xscore) ;
void printzq(double *z0q, double *z1q, double *z2q, int numqbins) ;
void ddadd(double ***ddcbins, double *z0q, double *z1q, double *z2q, int numqbins, int diffmax)  ; 
void ddcorr(CORR ***corrjbins, double ***ddcbins, int numdbins, int qbin) ;   
void fixjcorr(CORR **corrbins, CORR ***corrjbins, int numbins)  ;

void slocv(double *cout, double *a1, double *b1, int n, int maxlag) ;
void fftcv(double *cout, double *a1, double *b1, int n, int maxlag) ;
void fftauto(double *cout, double *a1, int n, int maxlag) ;
void fftcv2(double *cout, double *cout2, double *a1, double *b1, int n, int maxlag) ;
void dorunfit() ;
int dowork(char *admixpop, char **eglist) ; 

double dd0sum = 0 ; 
double dd0count = 0 ; 
int numdbins ; 

int main(int argc, char **argv)
{

  double ymem ;  
  int ret, i, copyret ; 
  int t, chrom ; 
  int nignore, numrisks = 1 ;
  char **eglist ;
  int numeg ; 
  char ***xjobs ; 
  int  njobs ; 
  char *dir ;
  char cwd[1024]  ;
  char sss[1024]  ;
  SNP *cupt ; 
  


  ofile = stdout; 
  packmode = YES ;
  

  printnl() ;
  printf("## DATES.  Version %s\n", WVERSION) ;
  printnl() ;

  readcommands(argc, argv) ;

  if (parname == NULL) return 0 ;

  printnl() ;
  cputime(0) ;
  calcmem(0) ;

  if (seed == 0) seed = seednum() ;
  SRAND(seed) ;
  printf("seed: %d\n", seed) ;

  if (snpname[0] != '/') { 
   ZALLOC(snppath, 256, char) ; 
   if (realpath(snpname, snppath) == NULL) fatalx("realpath can't find snpname!\n") ;
   printf("snppath: %s\n", snppath) ; 
   fflush(stdout) ; 
  }

  numsnps = 
    getsnps(snpname, &snpmarkers, 0.0, badsnpname, &nignore, numrisks) ;

  numindivs = getindivs(indivname, &indivmarkers) ;


  printf("genotypename:  %s\n", genotypename) ;
  fflush(stdout) ;

   getgenos(genotypename, snpmarkers, indivmarkers, 
    numsnps, numindivs, nignore) ;


   t = 0 ;
   for (i=0; i<numsnps; ++i) { 
    cupt = snpmarkers[i] ;
    cupt -> tagnumber = i ;
    chrom = cupt -> chrom ; 
    
    if (chrom<1) cupt -> ignore = YES ;
    if (chrom>22) cupt -> ignore = YES ;
    if ((xchrom>0) && (chrom != xchrom)) cupt -> ignore = YES ;
    if (chrom == xnochrom) cupt -> ignore = YES ;
    if (cupt -> ignore == NO) ++t ;
   }
   printf("valid snps: %d\n", t) ;

   numsnps = rmsnps(snpmarkers, numsnps, NULL) ;   
/**
   numindivs = rmindivs(snpmarkers, numsnps, indivmarkers, numindivs) ;
*/  

   if (checkmap && (cmap(snpmarkers, numsnps) == NO)) { 
    printf("running rolloff without a real map.  If you mean to do this set checkmap: NO\n") ;
    fatalx("no real map\n") ;
   }


  if (admixlist == NULL) { 
   ZALLOC(eglist, numindivs, char *) ; 
   numeg = loadlist(eglist, poplistname) ;
   if (numeg != 2) fatalx("must have exactly 2 admixing populations! poplistname: %s\n", poplistname) ;

   fflush(stdout) ;
   ret = dowork(admixpop, eglist) ; 

 }

  if (admixlist != NULL)  { 
   ZALLOC(eglist, 2, char *) ; 
   numeg = 2 ; 

   getcwd(cwd, 1024) ; 


   njobs = numlines(admixlist) ; 
   ZALLOC(xjobs, 4, char **)  ; 
   for (i=0; i<4; i++) { 
    ZALLOC(xjobs[i], njobs, char *) ; 
   } 

   njobs = getnames(&xjobs, njobs, 4, admixlist) ; 

   for (i=0; i<njobs; i++) {
    eglist[0] = xjobs[0][i] ;  
    eglist[1] = xjobs[1][i] ;  
    numeg = 2 ; 
    admixpop  = xjobs[2][i] ; 

    freestring(&oname) ; 
    freestring(&logfit) ; 
    
    dir = xjobs[3][i] ; 

//  mkdir(dir, 0777) ; 

    sprintf(sss, "mkdir -p -m 0777 %s\n", dir) ;
    system(sss) ; 

    sprintf(sss, "cp %s %s", parname, dir) ; 
    copyret = system(sss) ; printf("copy param code: %d\n", copyret) ;

    chdir (dir) ; 
    printf("calling dowork:  %s %s %s\n", eglist[0], eglist[1], admixpop) ;
    fflush(stdout) ;
    ret = dowork(admixpop, eglist) ;  
    chdir(cwd) ; 
   }
  }

  ymem = calcmem(1)/1.0e6 ;
  printf("##end of dates: %12.3f seconds cpu %12.3f Mbytes in use\n", cputime(1), ymem) ;
  return ret ;



}
int dowork(char *admixpop, char **eglist)  
{

  char s1[MAXSTR], s2[MAXSTR] ;
  int x, x2, i, j, k, g, s, t, d, tt, iter ; 
  int a, b, c, j1, j2, k1, k2  ;
  SNP *cupt, *cupt2 ;
  Indiv *indx, **indadmix ;

  double y, y1, y2, yy, ywt, z1, z2, znorm, yn, ytime ;
  double *yy2 ;
  int ind, chrom, lastchrom, printit, num0=0, numeg=2 ;
  double xc[9], xww[3] ;
  char sss[256] ;

  double *ttt ;
  double yscore, ybase ;
  double lambdabase ;
  double *scarray, *wacc, ywynn ;
  int ncount = 0 ;
  double lastscore = 0.0 ;
  double *tvecs, *pp ;
  double *w0, *w1, *w2, *res, alpha, *ww, *ww1, *ww2, *dd, *wt, *wk  ;
  double *wa, *wb ;  
  int *sindex, *xindex ; 
  double *xnum, *xden, *xdenq ;
  double *znum, *zden ;
  double **xxnum, **xxden ;
  int numbins  ;  
  double co1, co2, dis, yran ;     
  double wtmean, wtnum ; 
  CORR **corrbins ;
  CORR *corrpt  ;
  CORR ***corrjbins ;  //jackknife
  int nloop = 0 ;
  int nncount = 0 ;
  int nsnp1 = 0, nsnp2 = 0 ;
  int numchrom = 22 ;
  double **jwt ;
  int numx ;  // snps in main double loop
  double timeoffset ;  
 
  int numqbins ; 
  double *z0q, *z1q, *z2q ; 
  double **ddbins, ***ddcbins ; 
  int diffmax ; 
  double ymem ;  
  int xcnt[2] ;  // admixing pop sizes

  if (weightname != NULL) {
   printf("weights set on input\n") ;
  }
  if (runmode == 2) ldmode = YES ;
  if (zdipcorrmode) setzdipmode(1) ;

  if (chithresh < -0.5) {
   chithresh = 6.0 ;
  }


  printf("chithresh: %9.3f\n", chithresh) ;

  ZALLOC(admixpoplist, MAXAPOPS, char *) ;  
   if (admixpop == NULL) fatalx("no admixpop\n") ; 
    nadmixpop = 1 ;
    admixpoplist[0] = strdup(admixpop) ;

  if (oname == NULL) { 
    sprintf(sss, "%s.out", admixpop) ;  
    oname = strdup(sss) ;  
    printf("output tables: %s\n", oname) ;
  }


//   numsnps = rmsnps(snpmarkers, numsnps, NULL) ;   
  for (i=0; i<numsnps; i++) { 
   cupt = snpmarkers[i] ; 
   cupt -> ignore = NO ;
  }

   setfc(snpmarkers, numsnps) ; // does allocation
   numqbins = setqbins(snpmarkers, numsnps, binsize, qbin) ;       
   if (numqbins == 0) printf("no mesh\n") ; 
   else { 
    ZALLOC(qindex, numqbins, int) ;
    printf("mesh used :: qbin: %d  numqbins: %d\n", qbin, numqbins) ;
   }

   ncount = 0 ;
   ivzero(xcnt, 2) ;
   ZALLOC(indadmix, numindivs, Indiv *) ; 
   for (i=0; i<numindivs; ++i) { 
    indx = indivmarkers[i] ;
    indx -> affstatus = -1 ;
    indx -> ignore = NO ; 
    t = indxindex(eglist, numeg, indx -> egroup) ;
    if (t>=0) { 
     indx -> affstatus = t ;
     ++xcnt[t] ;  
    }
    tt = indxindex(admixpoplist, nadmixpop, indx -> egroup) ;
    if (tt>=0) { 
     if (t>0) fatalx(" %s in admixing and admixed pops\n", indx -> egroup) ;
     indx -> affstatus = 99 ;
     indadmix[ncount] = indx ;
     ++ncount ;   
    }
    if (indx -> affstatus  <  0) indx -> ignore = YES ;
    if (indx -> ignore) continue ;

   }
   
   if (ncount == 0) fatalx("no admixed samples found!\n") ;
   printf("number admixed: %d numeg: %d\n", ncount, numeg) ;
   for (k=0; k<numeg; ++k) { 
    printf("group %d %s samplesize: %4d\n", k, eglist[k], xcnt[k]) ;
   }
   if ((xcnt[0] == 0) || (xcnt[1] == 0)) fatalx("admixing population has no samples!\n") ;

//   numindivs = rmindivs(snpmarkers, numsnps, indivmarkers, numindivs) ;

   printf("numsnps: %d  numindivs: %d\n", numsnps, numindivs) ;
 
   if (numsnps == 0)   { 
     printf("zero snps -- misspelled population? -- exiting\n") ; 
     return -1 ; 
   }

   ZALLOC(w0, numsnps, double) ;
   ZALLOC(w1, numsnps, double) ;
   ZALLOC(w2, numsnps, double) ;
   ZALLOC(wk, numsnps, double) ;
   ZALLOC(ww1, numsnps  + numbins, double) ;
   ZALLOC(ww2, numsnps  + numbins, double) ;
   ZALLOC(wa, numsnps, double) ;
   ZALLOC(wb, numsnps, double) ;
   ZALLOC(res, numsnps, double) ;
   ZALLOC(dd, numsnps, double) ;
   ZALLOC(wt, numsnps, double) ;

   ZALLOC(sindex, numsnps, int) ;
   ZALLOC(xindex, numsnps, int) ;

   numbins = nnint(maxdis/binsize) + 5 ;
   ZALLOC(ww, numsnps +numbins, double) ;

   ZALLOC(xnum, numbins, double) ;
   ZALLOC(xden, numbins, double) ;
   ZALLOC(xdenq, numbins, double) ;
   ZALLOC(znum, numbins, double) ;
   ZALLOC(zden, numbins, double) ;
   ZALLOC(yy2, numbins, double) ;

    xxnum = initarray_2Ddouble(numchrom+1, numbins, 0) ;
    xxden = initarray_2Ddouble(numchrom+1, numbins, 1.0e-10) ;

    if (numqbins>0) { 
     ZALLOC(z0q, numqbins, double) ;  
     ZALLOC(z1q, numqbins, double) ;  
     ZALLOC(z2q, numqbins, double) ;  
     numdbins = numbins * qbin ; 
     ddbins = initarray_2Ddouble(7, numdbins, 0) ;   
     ZALLOC(ddcbins, numchrom+1, double **) ; 
     for (k=1; k<=numchrom; ++k) { 
      ddcbins[k] = initarray_2Ddouble(7, numdbins, 0) ;   
     }
    }

   ZALLOC(firstsnp, numchrom+1, int) ;
   ZALLOC(lastsnp, numchrom+1, int) ;
   ZALLOC(firstqb, numchrom+1, int) ; // first qbin in chromosome
   ZALLOC(lastqb, numchrom+1, int) ;

   vclear(ww, 1.0, numsnps) ;

   ZALLOC(corrbins, numbins, CORR *) ;

   for (x=0; x<numbins; ++x) { 
     ZALLOC(corrbins[x], 1, CORR) ;
     corrpt = corrbins[x] ; 
     clearcorr(corrpt) ;
   }

   ZALLOC(corrjbins, numchrom+1, CORR **) ;
   jwt = initarray_2Ddouble(numchrom+1, numbins, 0.0) ; 
   for (k=1; k<=numchrom; ++k) {
    ZALLOC(corrjbins[k], numbins, CORR *) ;

    for (x=0; x<numbins; ++x) { 
     ZALLOC(corrjbins[k][x], 1, CORR) ;
     corrpt = corrjbins[k][x] ; 
     clearcorr(corrpt) ;
    }
   }

   yn = (double) numsnps ;

/** 
 w1 is (normalized) parental allele freq difference 
 w2 is (normalized) (geno - expected)   where expected is freq average
*/

   for (x=0; x<numsnps; ++x)  { 
     
    cupt = snpmarkers[x] ;
    if (cupt -> ignore) continue ; 

    cntit1(xc, cupt, indivmarkers, numindivs, 0) ;  
     y1 = 2*xc[2] + xc[1] ;  
     y2 = 2*asum(xc,3) ;    
     if (y2<=0.01) { 
      cupt -> ignore = YES ;
      continue ;
     }
     wa[x] = cupt -> cauc_freq = y1/y2 ;

     cntit1(xc, cupt, indivmarkers, numindivs, 1) ;  
     y1 = 2*xc[2] + xc[1] ;  
     y2 = 2*asum(xc,3) ;    
     if (y2<=0.01) { 
      cupt -> ignore = YES ;
      continue ;
     }
     ++nsnp1 ; 
     wb[x] = cupt -> af_freq = y1/y2 ;
     y1 = cupt -> cauc_freq - cupt -> af_freq ;
     y = fabs(y1) ;    
     if (y<0.001) { 
      cupt -> ignore = YES ;
      continue ;
     }
     
     y = cupt -> cauc_freq + cupt -> af_freq ; 
     y *= 0.5 ;  
     y2 = y*(1.0-y) ; 
     if (y2<.001) {  
      cupt -> ignore = YES ;
      continue ;
     }
     ++nsnp2 ; 
     znorm = sqrt(y) ; // ?? involved in weight ?
     switch (runmode) { 
       case 0:  cupt -> weight = y1 / znorm ;
        break ;
       case 1:  cupt -> weight = y1 ; 
        break ; 
       case 2:  cupt -> weight = 1.0 ;
        break ;
       default: fatalx("bad runmode\n") ;
     }
     if (isnan(cupt -> weight)) fatalx("bad weight\n") ;
   }


   if (weightname != NULL) { 
    printf("loading weights...\n") ; 
    getweights(weightname, snpmarkers, numsnps) ;
   }

   if (timeoffsetname != NULL) { 
    printf("loading time offset...\n") ; 
    getindvals(timeoffsetname, indivmarkers, numindivs) ;

    ZALLOC(ttt, numindivs, double) ;
    for (k=0; k<ncount; ++k) { 
     indx = indadmix[k] ; 
     ttt[k] = timeoffset = indx -> qval ; 
    }
    y = asum(ttt, ncount) / (double) ncount ; 
    
    
    for (k=0; k<ncount; ++k) { 
     indx = indadmix[k] ; 
     timeoffset = indx -> qval ; 
     printf("timeoffset: %20s %9.3f\n", indx -> ID,  timeoffset) ;   
     indx -> qval -= y ; 
    } 
   }

   for (i=0; i < ncount; ++i) { 
// main loop 
    indx = indadmix[i] ;  
    timeoffset = indx -> qval ;
    k = indx -> idnum ; 
    x = 0 ; 
    ivclear(sindex, -1, numsnps) ;
    ivclear(firstsnp, 1000*1000*1000, numchrom+1) ;
    ivclear(lastsnp, -1, numchrom+1) ;
    ivclear(firstqb, 1000*1000*1000, numchrom+1) ;
    ivclear(lastqb, -1, numchrom+1) ;

    for (j=0; j<numsnps; ++j) {  
     cupt = snpmarkers[j] ;
     if (cupt -> ignore) continue ;  
     g = getgtypes(cupt, k) ;
     if (g<0) continue ;
     w0[x] = (double) g / 2.0 ; 
     w1[x] = cupt -> cauc_freq ;  
     w2[x] = cupt -> af_freq ;  
     wt[x] = w1[x] - w2[x] ;  // weight = freq diff ;  

     t = cupt -> chrom ; 

     firstsnp[t] = MIN(firstsnp[t], j) ; 
     lastsnp[t] =  MAX(lastsnp[t], j) ; 



     if (weightname != NULL) wt[x] = cupt -> weight ;  
     if (isnan(wt[x])) fatalx("bad weight!\n") ;

     sindex[j] = x ;
     xindex[x] = j ; 
     ++x ; 
    }
    vvm(ww1, w0, w2, x) ;
    vvm(ww2, w1, w2, x) ;

    y = vdot(ww1, ww2, x) / vdot(ww2, ww2, x) ;
    printf("coeff: %3d %20s %9.3f\n", i, indx -> ID, y) ;
    vst(ww1, w1, y, x) ; 
    vst(ww2, w2, 1.0-y, x) ; 
    vvp(ww1, ww1, ww2, x) ; // prediction
    vvm(res, w0, ww1, x) ;  // residual 
    numx = x ; 

    if (numx==0) fatalx("no valid snps!\n") ;


    vzero(z0q, numqbins) ;
    vzero(z1q, numqbins) ;
    vzero(z2q, numqbins) ;

   if (numqbins>0)  {

     ivclear(qindex, -1, numqbins) ;

     for (k1=0; k1<numx; ++k1)  { 

      j1 = xindex[k1] ;
      cupt = snpmarkers[j1] ;

      if (cupt -> ignore) continue ;
      s = cupt -> tagnumber  ; 
      y = res[k1] * wt[k1]  ; 
      qindex[s] = k1 ; // just 1 snp maps

      t = cupt -> chrom ; 

      firstqb[t] = MIN(firstqb[t], s) ;
      lastqb[t] =  MAX(lastqb[t], s) ;

      z0q[s] += 1 ; 
      z1q[s] += y ; 
      z2q[s] += y*y ; 
    } 

   y = (double) qbin * maxdis / binsize ;
   diffmax = nnint(y) ; 
   ddadd(ddcbins, z0q, z1q, z2q, numqbins, diffmax) ;

  }

  if (numqbins == 0) printf("starting main loop. ID %s numsnps: %d\n", indx -> ID, numsnps) ;
  lastchrom = -1 ; 
  for (k1=0; k1<numx; ++k1)  { 
   if (numqbins > 0) break ; 
   j1 = xindex[k1] ;
   cupt = snpmarkers[j1] ;
   if (cupt -> ignore) continue ;
   chrom = cupt -> chrom ;
   printit = -1 ;
   if (lastchrom == -1) lastchrom = chrom ; 
   if (lastchrom != chrom) { 
    printit = chrom ; 
    lastchrom = chrom ;
   }
    for (k2=k1+1; k2 < numx; ++k2) { 
     j2 = xindex[k2] ;
     cupt2 = snpmarkers[j2] ;
     if (cupt2 -> ignore) continue ;
      t = cupt2 -> chrom - cupt -> chrom ;
      if (t != 0) break ;
      dis = cupt2 -> genpos - cupt -> genpos ;
      if (dis<0.0) fatalx("badbug\n") ;
      if (dis >= maxdis) break ;
      s = (int) (dis/binsize) ; // as in simpsim2 
      ytime = exp(-timeoffset*dis) ;

      yy  = res[k1]*res[k2] ;                    
      ywt = wt[k1]*wt[k2] ;                     

      znum[s] += ywt * yy ;  
      zden[s] += ywt * ywt ;

      y = exp(-timeoffset*dis) ;

      ywt *= y ;

      xnum[s] += ywt * yy ;  
      xden[s] += ywt * ywt ;

      t = ranmod(1000*1000) ; 
      if (t==-1)  {  
       printf("zza %12.6f ", dis) ;  
       printf("%d %12.6f ", s, y) ;  
       y1 = (double) s * binsize ;  
       printf("%12.6f ", exp(-timeoffset*y1) ) ; 
       printf(" %12.6f %12.6f", xnum[s]/znum[s]) ;  
       printf(" %12.6f %12.6f", xden[s]/zden[s]) ;  
       printnl() ;
      }
  

      yy2[s] += yy * yy ;
      y1 = res[k1]*wt[k1] ;
      y2 = res[k2]*wt[k2] ;
      y1 *= ytime ;
      
      corrpt = corrbins[s] ; 
      addcorr(corrpt, y1, y2) ;
      ++jwt[chrom][s] ;
      for (k=1; k<= numchrom; ++k) {
        if (k==chrom) continue ;
        if (jackknife == NO) break ;
        corrpt = corrjbins[k][s] ;
        addcorr(corrpt, y1, y2) ;
        xxnum[k][s] += ywt*yy ; 
        xxden[k][s] += ywt*ywt ; 
      }
   }
   if (printit > 0) { 
    printf("sample: %20s chrom: %3d done\n", indx -> ID, printit-1) ;  fflush(stdout) ;
   }
  }

  
   
   if (numqbins==0) printf("sample: %20s chrom: %3d done\n", indx -> ID, lastchrom) ;  fflush(stdout) ;
   y = cputime(1) ;  
   printf("sample: %20s completed.  Total CPU since start: %9.3f\n", indx -> ID, y) ; fflush(stdout) ; 
 }

 if (numqbins>0) {  

   ddcorr(corrjbins, ddcbins,  numdbins, qbin) ;  
   fixjcorr(corrbins, corrjbins, numbins) ; // set up jackknifed corrs 

   vzero(ww, numbins) ;
   sprintf(sss, " ##Z-score and correlation:: %s  binsize: %12.6f", admixpop, binsize) ;
   dumpit(oname, ww, corrbins, numbins-5, binsize, sss) ;

   for (k=1; k<=numchrom; ++k) { 
    if (jackknife == NO) break ;
    sprintf(s1, "%s:%d", oname, k) ;
    sprintf(s2, "## Jackknife output: chrom %d", k) ;
    dumpit(s1, ww, corrjbins[k], numbins-5, binsize, s2) ;
   }

 } 


  if (numqbins == 0) {

  vsp(xden, xden, 1.0e-20, numbins) ; 
  vsp(zden, zden, 1.0e-20, numbins) ; 
  vsp(yy2, yy2, 1.0e-10, numbins) ; 

   vsqrt(xden, xden, numbins) ;

   vvd(ww, xnum, xden, numbins) ;
   vvd(ww1, znum, zden, numbins) ;
   
   sprintf(sss, " ##Z-score and correlation:: %s  binsize: %12.6f", poplistname, binsize) ;
   dumpit(oname, ww, corrbins, numbins-5, binsize, sss) ;
  
   for (k=1; k<=numchrom; ++k) { 
    if (jackknife == NO) break ;
    sprintf(s1, "%s:%d", oname, k) ;
    sprintf(s2, "## Jackknife output: chrom %d", k) ;
    vvd(ww, xxnum[k], xxden[k], numbins) ;
    dumpit(s1, ww, corrjbins[k], numbins-5, binsize, s2) ;
   }
  }

  dorunfit() ;

   free(w0) ;
   free(w1) ;
   free(w2) ;
   free(wk) ;
   free(ww1) ;
   free(ww2) ;
   free(wa) ;
   free(wb) ;
   free(res) ;
   free(dd) ;
   free(wt) ;

   free(sindex) ;
   free(xindex) ;

   free(ww) ;

   free(xnum) ;
   free(xden) ;
   free(xdenq) ;
   free(znum) ;
   free(zden) ;
   free(yy2) ;

    free2D(&xxnum, numchrom+1) ;
    free2D(&xxden, numchrom+1) ;

    if (numqbins>0) { 
     free(z0q) ; 
     free(z1q) ; 
     free(z2q) ; 

     for (k=1; k<=numchrom; ++k) { 
      free2D(&ddcbins[k], 7) ;
     }
     free(ddcbins) ; 
     free2D(&ddbins, 7) ; 
    }

   free(firstsnp) ; 
   free(lastsnp) ; 
   free(firstqb) ; 
   free(lastqb) ; 


   for (x=0; x<numbins; ++x) { 
     free(corrbins[x]) ;
   }

   free(corrbins) ;


   free2D(&jwt, numchrom+1) ;

   for (k=1; k<=numchrom; ++k) {

    for (x=0; x<numbins; ++x) { 
     free(corrjbins[k][x]) ; 
    }
    free(corrjbins[k]) ; 
   }

   free(corrjbins) ;


  return 0;

}

void dorunfit() 
{  
  char s1[256], s2[256], s3[256], sss[256] ;

  if (runfit == NO) return ; 
   if (logfit == NULL) { 
    sprintf(s2, "%s:log", admixpop) ; 
    logfit = strdup(s2) ; 
  }

   if (jackknife) {
     sprintf(s1, "dates_jackknife -p %s -l %9.3f -z %s -r %d ", parname, lovalfit, admixpop, seed) ;
     if (snppath != NULL) { 
      sprintf(s3, "-m %s", snppath) ; 
      strcat(s1, s3) ;
     }
   }
   else  {
     sprintf(s1, "run_dates_expfit -p %s -l %9.3f -z %s -r %d ", parname, lovalfit, admixpop, seed) ;
   }
   
   if (afffit) { 
    strcat(s1, " -a") ; 
   }

   sprintf(s2, " >%s", logfit) ;
   strcpy(sss, s1) ; 
   strcat(sss, s2) ; 
   printf("calling fit!...\n") ; 
   printf("%s\n", sss) ; 
   fflush(stdout) ;
   system(sss) ;

}

void dumpit(char *oname, double *ww, CORR **corrbins, int len, double bsize, char *hdr) 
{
   double y, yreg, ycovar ; 
   int k, ret ;
   FILE *fff ;
   CORR *corrpt ;

   if (oname == NULL) return ;
   openit(oname, &fff, "w") ;

   if (hdr != NULL) fprintf(fff, "%s\n", hdr) ;
   for (k=0; k<len; ++k) { 

    y = (k+1) * bsize ;
    corrpt = corrbins[k] ;
    fprintf(fff, "%9.3f ", 100.0*y) ;  // CM
//  fprintf(fff, "%12.6f ", ww[k]) ;   
//  fprintf(fff, "%12.6f ", corrpt -> Z) ;
    if (runmode != 2) {
     calccorr(corrpt, 1, NO) ;
     ycovar = corrpt -> v12 ; 
     yreg = corrpt -> S12 / (corrpt -> S11 + 1.0e-20) ;
     if (isnan(yreg)) { 
      printf("bad corrpt!\n") ; 
      printcorr(corrpt) ; 
      fatalx("dumpit fails\n") ; 
     }
     fprintf(fff, "%15.9f ", ycovar) ; 
     fprintf(fff, "%12.6f ", yreg) ; 
     fprintf(fff, "%12.6f ", corrpt -> corr) ;   
     fprintf(fff, "%12.0f ", corrpt -> S0) ;
     if (corrpt -> S0 > 1.0e20) { 
      printf("bad corr!\n") ; 
      printcorr(corrpt) ; 
      fatalx("bad corr!\n") ; 
     }
     if (fabs(yreg) > 10.0)  { 
      printf("bad corr!\n") ; 
      printcorr(corrpt) ; 
      fatalx("bad yreg!\n") ; 
     }
    }
    fprintf(fff, "\n") ;

   }


   fclose(fff) ;

}
void
dorcx(double *ans, double *rr, double *dd, double binsize, double maxdis) 
// accumulate dot product into bins
{
   int x, x2, z, t, s, k ; 
   SNP *cupt, *cupt2 ;
   double dis ;
   double *xnum, *xd1, *xd2, y1, y2 ;
   double *rrr ;

   z =nnint(maxdis/binsize) ;  
   ZALLOC(xnum, z+1, double) ;
   ZALLOC(xd1, z+1, double) ;
   ZALLOC(xd2, z+1, double) ;

   vclear(xd1, 1.0e-10, z+1) ;
   vclear(xd2, 1.0e-10, z+1) ;

   for (x=0; x<numsnps; ++x)  {
    cupt = snpmarkers[x] ; 
    for (x2=x+1; x2<numsnps; ++x2)  {
      cupt2 = snpmarkers[x2] ; 
      t = cupt2 -> chrom - cupt -> chrom ;
      if (t != 0) break ;
      dis = cupt2 -> genpos - cupt -> genpos ;
      if (dis >= maxdis) break ;
      s = nnint(dis/binsize) ;
      y2 = dd[x]*dd[x2] ;
      for (k=0; k<numindivs; ++k) { 
       rrr = rr + k*numsnps ;

       switch (runmode) { 
        case 0:  
         y2 = dd[x]*dd[x2] ;
         y1 = rrr[x]*rrr[x2] ;
         break ;
      
        case 1:  
         y2 = rrr[x2]*dd[x2] ;
         y1 = rrr[x]*dd[x]  ;
         break ;

/**
        case 2:  
         y2  = 1.0 ;
         y1 = rrr[x]*rrr[x2]  ;
         break ;

        case 3:  
         y2  = 1.0 ;
         y1 = dd[x]*dd[x2] ;
         break ;
*/

        default: 
         fatalx("bad runmode\n") ; 
       }

       if (isnan(y1)) fatalx("bad y1\n") ;
       if (isnan(y2)) fatalx("bad y2\n") ;
       xnum[s] += y1*y2 ;
       xd1[s] += y1*y1 ;
       xd2[s] += y2*y2 ;
      }
    }
   }
  
   
   for (s=0; s <= z; ++s) {
    ans[s] = xnum[s]/sqrt(xd1[s]*xd2[s]) ;
   }

   free(xnum) ;
   free(xd1) ;
   free(xd2) ;
}

void
dorc2(double *ans, double *res, double binsize, double maxdis) 
// accumulate square of dot product into bins
{
   int x, x2, z, t ; 
   SNP *cupt, *cupt2 ;
   double dis, y ;

   z =nnint(maxdis/binsize) ;  
   vzero(ans, z+1) ;
   for (x=0; x<numsnps; ++x)  {
    cupt = snpmarkers[x] ; 
    for (x2=x+1; x2<numsnps; ++x2)  {
      cupt2 = snpmarkers[x2] ; 
      t = cupt2 -> chrom - cupt -> chrom ;
      if (t != 0) break ;
      dis = cupt2 -> genpos - cupt -> genpos ;
      if (dis >= maxdis) break ;
      z = nnint(dis/binsize) ;
      y  = res[x]*res[x2] ;
      ans[z] += y*y ;             
    }
   }
}

int readcommands(int argc, char **argv) 

{
  int i,haploid=0;
  phandle *ph ;
  char str[5000]  ;
  char *tempname ;
  int n ;

  while ((i = getopt (argc, argv, "p:vV")) != -1) {

    switch (i)
      {

      case 'p':
	parname = strdup(optarg) ;
	break;

      case 'v':
	printf("version: %s\n", WVERSION) ; 
	break; 

      case 'V':
	verbose = YES ;
	break; 

      case '?':
	printf ("Usage: bad params.... \n") ;
	fatalx("bad params\n") ;
      }
  }

         
   if (parname == NULL) return -1 ;
   printf("parameter file: %s\n", parname) ;
   ph = openpars(parname) ;
   dostrsub(ph) ;

   getstring(ph, "genotypename:", &genotypename) ;
   getstring(ph, "snpname:", &snpname) ;
   getstring(ph, "indivname:", &indivname) ;
   getstring(ph, "poplistname:", &poplistname) ;
   getstring(ph, "weightname:", &weightname) ;
   getstring(ph, "timeoffset:", &timeoffsetname) ;
   getstring(ph, "admixpop:", &admixpop) ;
   getstring(ph, "admixlist:", &admixlist) ;
   getstring(ph, "badsnpname:", &badsnpname) ;
   getstring(ph, "dumpname:", &dumpname) ;
   getstring(ph, "loadname:", &loadname) ;
   getint(ph, "popestimate:", &popestimate) ;
   getint(ph, "ransample:", &ransample) ;
   getint(ph, "chrom:", &xchrom) ;
   getint(ph, "nochrom:", &xnochrom) ;
   getstring(ph, "output:",   &oname) ;
   getint(ph, "runmode:", &runmode) ;
   getint(ph, "ldmode:", &ldmode) ;
   getint(ph, "seed:", &seed) ;
   getint(ph, "nxlim:", &nxlim) ;
   getint(ph, "minparentcount:", &minparentcount) ;
   getint(ph, "norun:", &norun) ;
   getint(ph, "flatmode:", &flatmode) ;
   getint(ph, "zdipcorrmode:", &zdipcorrmode) ;
   getint(ph, "jackknife:", &jackknife) ;

   getdbl(ph, "maxdis:", &maxdis) ;
   getdbl(ph, "binsize:", &binsize) ;
   getdbl(ph, "chithresh:", &chithresh) ;
   getint(ph, "checkmap:", &checkmap) ;
   getint(ph, "qbin:", &qbin) ;
// parameters for runfit etc. 
   getint(ph, "runfit:", &runfit) ;
   getint(ph, "afffit:", &afffit) ;
   getdbl(ph, "lovalfit:", &lovalfit) ;

   writepars(ph) ;
   closepars(ph) ;

}

int      
setqbins(SNP **snpm, int numsnps, double binsize, int qbin) 
// return number of qbins ;  set tagnumber to qbin 
{
 double qb, ydis, lastgenpos ; 
 int chrom, s = 0, k ; 
 SNP *cupt ;

 if (qbin==0) return 0 ; 

 qb = binsize/(double) qbin ; 

 for (k=0; k<numsnps; ++k) { 
  cupt = snpm[k] ; 
  if (k==0) { 
   ydis = 0 ; 
   lastgenpos = cupt -> genpos ;
   s = floor(ydis/qb) ;  
   cupt -> tagnumber = s ; 
   chrom = cupt -> chrom ; 
   continue ; 
  }
  if (cupt -> chrom != chrom) { 
   ydis += 5 ; 
   chrom = cupt -> chrom ; 
  }
  else { 
   ydis += (cupt -> genpos - lastgenpos) ;
  }
  lastgenpos = cupt -> genpos ;
  s = floor(ydis/qb) ;  
  cupt -> tagnumber = s ;
 }

 return s+1 ;



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
    if (cupt -> ignore) continue ;
    chrom = cupt -> chrom ;
//  if ((noxdata) && (chrom==23)) continue  ;
    fchrom[chrom] = MIN(fchrom[chrom],i) ;
    lchrom[chrom] = MAX(lchrom[chrom],i) ;
  }
}
double cntit1(double *xc, SNP *cupt, Indiv **indm, int numindivs, int t) 
{
  Indiv *indx ;
  int k, g ;

  vzero(xc, 3) ;
  if (cupt -> ignore) return -1 ;
  for (k=0; k<numindivs; ++k) { 
   indx = indm [k] ;
   if (indx -> ignore) continue ;
   if (indx -> affstatus  != t) continue ;
   g = getgtypes(cupt, k) ;
   if (g<0) continue ;
   ++xc[g] ;
  }
  return asum(xc,3) ;
}


double cntit2(double *xc, SNP *cupt, SNP *cupt2, Indiv **indm, int numindivs, int t) 
{
  Indiv *indx ;
  int k, e, f ;

  vzero(xc, 9) ;
  if (cupt -> ignore) return -1 ;
  if (cupt2 -> ignore) return -1 ;
  for (k=0; k<numindivs; ++k) { 
   indx = indm [k] ;
   if (indx -> ignore) continue ;
   if (indx -> affstatus  != t) continue ;
   e = getgtypes(cupt, k) ;
   if (e<0) continue ;
   f = getgtypes(cupt2, k) ;
   if (f<0) continue ;
   ++xc[3*e+f] ;
  }
  return asum(xc,9) ;

}

void cx1(double *xc, double *y1, double *y2) 
{
   double x1[3], x2[3] ;
   int e, f ;

   vzero(x1,3) ;
   vzero(x2,3) ;
   for (e=0; e<3; ++e) {  
    for (f=0; f<3; ++f) {  
     x1[e] += xc[3*e+f] ;
     x2[f] += xc[3*e+f] ;
    }
   }
  
   *y1 = cminor(x1) ;
   *y2 = cminor(x2) ;
    
}

double cminor(double *xx) 
{
  double y1, y2 ;
  y1 = xx[1]+2*xx[2] ; 
  y2 = 2*asum(xx, 3) - y1 ;
  return MIN(y1, y2) ;
}
void
printzq(double *z0q, double *z1q, double *z2q, int numqbins)
{
  int s, t ;

  printf("##printzq %9.0f\n", asum(z0q, numqbins) ) ;

  for (s=0; s<numqbins; ++s)  { 
   t = nnint(z0q[s]) ;  
   if (t==0) continue ;  
   printf("zq %9d %6.0f %9.3f %9.3f\n", s, z0q[s], z1q[s], z2q[s]) ;
 }

}


void
ddadd(double ***ddcbins, double *z0q, double *z1q, double *z2q, int numqbins, int diffmax) 

// this is the only potentially quadratic part of the code.   
{
  double *dd00, *dd01, *dd10, *dd11, *dd02, *dd20 ; 
  double **ddbins ;
  double a0, a1, a2, b0, b1, b2, y ; 
  int s, t, tmax, d, js, jt, k, slo, shi, len, j ; 
  SNP *cupts, *cuptt ;

  double *dans ;  
  double *dans2 ;  

  if (numqbins == 0) return ; 
  ZALLOC(dans, diffmax +1, double) ;
  ZALLOC(dans2, diffmax +1, double) ;
 

 for (k=1; k<=numchrom; ++k) { 

  slo = firstqb[k] ;
  shi = lastqb[k] ;

  ddbins = ddcbins[k] ;

  dd00 = ddbins[0] ; 
  dd01 = ddbins[1] ; 
  dd10 = ddbins[3] ; 
  dd11 = ddbins[4] ; 
  dd02 = ddbins[2] ; 
  dd20 = ddbins[6] ; 

  len = shi-slo +1 ;

  if (len <= 1) continue ;

  s = slo ; 
  
  fftauto(dans, z0q +s, len, diffmax) ; 
  vvp(dd00, dd00, dans, diffmax+1) ;  
  for (j=0; j<=diffmax; ++j) { 
    y = fabs(dans[j]) ;
    if (isnan(y)) fatalx("bad ddadd\n") ;
    if (y> 1.0e20) { 
     printmat(z0q+s, 1, len) ; 
     fatalx("bad ddadd\n") ;
    }
  }

  fftcv2(dans, dans2, z0q +s, z1q +s, len, diffmax) ; 

  vvp(dd01, dd01, dans, diffmax+1) ;  
  vvp(dd10, dd10, dans2, diffmax+1) ;  

  for (j=0; j<=diffmax; ++j) { 
    if (isnan(dans[j])) fatalx("bad ddadd\n") ;
    if (isnan(dans2[j])) fatalx("bad ddadd\n") ;
  }

  fftauto(dans, z1q +s, len, diffmax) ; 
  vvp(dd11, dd11, dans, diffmax+1) ;  

  fftcv2(dans, dans2, z0q +s, z2q +s, len, diffmax) ; 
  vvp(dd02, dd02, dans, diffmax+1) ;  
  vvp(dd20, dd20, dans2, diffmax+1) ;  


 }
 
 free(dans) ;
 free(dans2) ;

}

void ddcorr(CORR ***corrjbins, double ***ddcbins, int numdbins, int qbin)  
{
  double *dd00, *dd01, *dd10, *dd11, *dd02, *dd20 ; 
  int s, d, k ; 
  CORR *corrpt, **corrtemp ; 
  double ys, dbinsize ; 
  double **ddbins, **ddtemp ;   

  if (numdbins==0) return ; 

  ddbins = initarray_2Ddouble(7, numdbins, 0) ;

  for (k=1; k<=numchrom; ++k) { 
   corrtemp = corrjbins[k] ; 
   ddtemp = ddcbins[k] ; 

  dd00 = ddtemp[0] ; 
  dd01 = ddtemp[1] ; 
  dd10 = ddtemp[3] ; 
  dd11 = ddtemp[4] ; 
  dd02 = ddtemp[2] ; 
  dd20 = ddtemp[6] ; 

  dbinsize = binsize / (double) qbin ; 
  for (d=1; d<numdbins; ++d) {   // dangerous bend,  lag 0 should be 0 and is not computed correctly!
   if (dd00[d] < 0.5) continue ; 
   ys = (double) d  * dbinsize ;     
   s = (int) (ys/binsize) ; 
   corrpt = corrtemp[s] ; 
   addcorr2(corrpt, dd00[d], dd01[d], dd10[d], dd11[d], dd02[d], dd20[d]) ; 

  }
 
 }

}
void fixjcorr(CORR **corrbins, CORR ***corrjbins, int numbins) 
{
  int k, x ; 
  CORR *corrpt ; 
 
  for (x=0; x<numbins; ++x) { 
// step 1 calculate global corr
    corrpt = corrbins[x] ;  
    clearcorr(corrpt) ; // not really necessary 
    for (k=1; k<=numchrom; ++k) { 
     pluscorr(corrpt, corrpt, corrjbins[k][x]) ; 
    }
// step 2 calculate jackknifed correlations by subtraction
    for (k=1; k<=numchrom; ++k)  {  
     minuscorr(corrjbins[k][x], corrpt, corrjbins[k][x]) ;
    }
  }

}
