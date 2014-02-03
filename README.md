DAIRRy-BLUP
============

A Distributed AI-REML Ridge Regression Best Linear Unbiased Prediction framework for genomic prediction.

This software was developed by Arne De Coninck and can only be used for research purposes.

Genomic datasets used for genomic prediction are constantly growing due to the decreasing costs of genotyping and increasing interest in improving agronomic performance of animals and plants. To be able to deal with those large-scale datasets, a distributed-memory framework was developed based on a message passing interface and a ScaLAPACK library. The complexity of the algorithm is defined by the number of SNPs included in the genomic prediction setting; the number of individuals only has a linear effect on the read-in time. To enhance performance it is advised to compile and execute DAIRRy-BLUP on an MPI-optimized machine.

#Installation

## Dependencies

DAIRRy-BLUP relies heavily on the following software packages, which have to be installed prior to installation of DAIRRy-BLUP. These software packages are all open source, except for the vendor-optimized implementations.

1. MPI ([OpenMPI](http://www.open-mpi.org/), [MPICH](http://www.mpich.org/), [IntelMPI](http://software.intel.com/en-us/intel-mpi-library))
2. [ScaLAPACK](http://www.netlib.org/scalapack/) and all its dependencies BLAS, BLACS, LAPACK, PBLAS (It is recommended to install a [vendor optimized implementation](http://www.netlib.org/scalapack/faq.html#1.3) )
3. HDF5 (http://www.hdfgroup.org/HDF5/)
4. CMake (http://www.cmake.org/)

Currently, compilation will only work with the Intel MKL libraries installed. When MKL libraries are not available, one must change the MKL libraries in the CMakelists.txt file to the ones which are installed. 

## Step-by-step

1. Unpack zip-file or clone git-repository
2. go into the directory `DAIRRy-BLUP`
3. make a new directory `build`
4. go into the directory `build`
5. type `cmake ..`
6. type `make`

# Usage

DAIRRy-BLUP only needs an input file to start. A default input file is provided: `defaultinput.txt`, more information on the arguments in the input-file can be found on the wiki.

To launch DAIRRy-BLUP with the default input file using for example 4 processes, the following command should be entered in the `build` directory:
`mpirun -np 4 ./DAIRRy-BLUP defaultinput.txt`

# Output

DAIRRy-BLUP creates 1 or 2 output-files depending on whether a test-dataset is given in the input file.
* `estimates.txt`: Lists the estimates for the fixed and random effects. First t values are estimates for the t fixed effects. Last s values are estimates of the s SNP effects 
* `EBV.txt`: When a test dataset is provided, breeding values are estimated for the individuals based on the estimated fixed and SNP effects and th egenotypes of the test individuals.

# Version history

* Version 0.1 (12/2013):
  1. First public release of DAIRRy-BLUP

# Contact

Please feel free to contact arne.deconinck[at]ugent.be for any questions or suggestions. 

[![Analytics](https://ga-beacon.appspot.com/UA-29101865-2/DAIRRy-BLUP/page)](https://github.com/igrigorik/ga-beacon)
