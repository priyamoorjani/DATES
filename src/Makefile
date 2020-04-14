HOMEL=$(PWD)    
BIN=$(HOMEL)//bin
IDIR=-I.   -Inicksrc


ND = ./nicksrc
NLIB = $(ND)/libnick.a




override LDLIBS +=  -lm -lgsl -lblas -lfftw3 $(NLIB)  
override CFLAGS += -c -g -p -Wimplicit  $(IDIR)  

## Harvard O2 
override CFLAGS += -I/n/app/openblas/0.2.19/include -I/n/app/gsl/2.3/include -I/n/app/fftw/3.3.7/include
override LDFLAGS += -L/n/app/openblas/0.2.19/lib -L/n/app/gsl/2.3/lib/  -L/n/app/fftw/3.3.7/lib



override CFLAGS += -g -pg -DHAVE_CONFIG_H    -Wimplicit-int -L. -L  -g -pg -D_IOLIB=2


$(NLIB): 
	$(MAKE) -C $(ND) 

nicklib: $(NLIB) 
	ln $(NLIB) nicklib

all:   dates_expfit dates  simpjack2 dowtjack grabpars
dates_expfit:    nicklib dates_expfit.o   regsubs.o  fitexp.o   gslfit.o
dates: $(NLIB) dates.o qpsubs.o mcio.o ldsubs.o admutils.o egsubs.o regsubs.o fftsubs.o
dates: $(NLIB) dates.o qpsubs.o mcio.o ldsubs.o admutils.o egsubs.o regsubs.o fftsubs.o
dowtjack: $(NLIB) dowtjack.o
simpjack2: $(NLIB) simpjack2.o 

grabpars: $(NLIB) grabpars.o  

clean:  
	rm -f  *.o  nicksrc/*.o nicklib

clobber:  clean
	rm -f  dates_expfit dates simpjack2 dowtjack grabpars nicksrc/libnick.a bin/*

install: all
	cp dates_expfit dates grabpars simpjack2 dowtjack perlsrc/*  ./bin 
