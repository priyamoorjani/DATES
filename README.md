# DATES

DATES (Distribution of Ancestry Tracts of Evolutionary Signals) is a method to estimate the time of admixture in ancient DNA samples described in Narasimhan, Patterson et al. 2018

#### Installation
To make the DATES executable, you will need the following libraries:
```
GSL 2.5 (https://www.gnu.org/software/gsl/)
OpenBLAS 0.2.20 (www.openblas.net)
FFTW 3.3.8 (http://www.fftw.org/download.html)
````
You will need to edit the file ``Makefile`` and change the following lines to the appropriate directories containing the ``include`` and ``lib`` files once you have these libraries installed.
```
override CFLAGS += -I<PATH_TO_GSL_INCLUDE> -I<PATH_TO_OPENBLAS_INCLUDE> -I<PATH_FFTW_INCLUDE>
override LDFLAGS += -I<PATH_TO_GSL_LIB> -I<PATH_TO_OPENBLAS_LIB> -I<PATH_FFTW_LIB>
```

 Once the ``Makefile`` is updated, simply run ``make`` and then add the executables in the ``bin`` directory to your PATH. 

```
make install
export PATH=$PATH:<PATH_TO_bin_directory>
```
Note, the ``bin`` directory contains additional executables like ``dates_jackknife`` which are needed for the jackknife and so skipping this step will lead to errors.

#### Input
DATES requires that the input data is available in one of these formats (See https://reich.hms.harvard.edu/software/InputFileFormats). To convert to the appropriate format, one can use ``CONVERTF`` program (See https://github.com/argriffing/eigensoft/tree/master/CONVERTF for details). 

#### Command line 
```
./dates -p $parfile >$logfile
```
$logfile: Name of the logfile. The DATES program prints various statistics to standard output which should be directed to the logfile.  <br />
$parfile: Name of parameter file.  <br />

#### Parameter file
```
genotypename: input genotype filename   # in eigenstrat format
snpname:    input snp filename          # in eigenstrat format
indivname:  input indiv filename        # in eigenstrat format
admixlist:  filename                    # This file contains the source and admixed populations to use for the analysis. Each line has the format: <source1> <source2> <testpopulation> <output_directory>, where source1 and source2 are the reference populations for the ancestral populations, testpop is the name of the admixed population and the output_directory is the name of the output directory.Output files are of the format output_directory/testpopulation.out
binsize:    number                      # in Morgans, range is from 0-1. Optimal binsize of 0.001 is recommended.
maxdis:     number                      # in Morgans, range is 0-1. For quicker runs, use max_distance < 1.0. However, for recent admixture,   ensure that max_distance is greater than the expected admixture LD blocks.
seed:       number                      # Random seed to ensure reproducibility of runs. 
jackknife: YES/ NO                      # if YES, program will run jackknife by dropping one chr in each run and estimate the mean and standard error across the 23 runs.
qbin:       number                      # discretization parameter on mesh size for the binned residuals. Higher qbin correlates loosely with higher accuracy and highly with longer run time.
runfit:  YES/NO                         # run exponential fit using least squares on the output to infer the date of admixture?
afffit:    YES/NO                       # use affine for the fit? 
lovalfit:  0.45                         # in centiMorgans, starting genetic distance.
```

##### Optional paramaters
```
weightname: weight_file            # Contains a weight for each SNP to be included in the run. If this parameter is not specified, the program uses the allele frequency differentiation between the ancestral populations as the weight for each SNP. 
minparentcount: number             # Contains the minimum number of ancestral individuals from each population that is required for the run. Default = 10.
chrom:      chromosome_number      # The analysis is limited to the specified chromosome only.
nochrom:    chromosome_number      # The specified chromosome is excluded from the analysis.
badsnpname: badsnp_list            # File contains a list of SNPs to be excluded from the analysis. 
```

#### Output
The program generates several output files (and some tempfiles).
```
output.out                        # This file contains output for the entire genome. The estimates of covariance values at various genetic distances, binned according to input values.
output.out:$chr, where chr=1-22   # These files contain the output for the jackknife where we remove one chromosome ($chr) in each run.
output.jin                        # Summary of the jackknife output. Columns are: <chr> <# SNPs on the chromosome> <Estimate date by removing the chromosome>
output.jout                       # Final output. <mean time of admixture> <SE based on jackknife>
output.fit                        # output of least squares exponential fit. <genetic_distance_inCM> < output> < fitted_value> <output-fitted output>
output.pdf                        # pdf of output with exponential fit.
```
#### Example 
An example run is available in ``src/example/`` directory. The data was simulated using ancestral haplotypes from 1000 Genomes Project West Africans (YRI) and Northern Europeans (CEU) using the our ADMIX simulator (https://github.com/priyamoorjani/Admix_simulator). For details of the run, please see ``post.sh`` and parfiles (``par.simulation, par.dates``) in the ``example/`` directory. Logfiles are provided for reference.

#### Support
Send queries to Priya Moorjani (moorjani@berkeley.edu) or Nick Patterson (nickp@broadinstitute.org)
