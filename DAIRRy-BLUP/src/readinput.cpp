#include <stdio.h>
#include <stdlib.h>
#include <iosfwd>
#include <string>
#include "shared_var.h"
#include <shared_var.h>
#include <iostream>
#include <fstream>
using namespace std;

int ri(char * filename) {
    std::ifstream inputfile(filename);
    string line;
    bool obs_bool=false, SNP_bool=false, fixed_bool=false, fixedfile_bool=false, datafile_bool=false, lam_bool=false, phenopath_bool=false, testpath_bool=false;
    bool eps_bool=false, maxit_bool=false, blocksize_bool=false, testfile_bool=false, ntest_bool=false, h5_bool=false, genopath_bool=false;

    fDn=( char* ) calloc ( 100,sizeof ( char ) );
    fXn=( char* ) calloc ( 100,sizeof ( char ) );
    Ts=( char* ) calloc ( 100,sizeof ( char ) );
    Sd = ( char* ) calloc ( 100,sizeof ( char ) );
    pd = ( char* ) calloc ( 100,sizeof ( char ) );

    l=100;
    b=64;
    e=0.01;
    mi=20;
    nt=0;
    Ccop=0;


    while( std::getline (inputfile,line)) {
        if(line=="#Observations") {
            std::getline (inputfile,line);
            n=atoi(line.c_str());
            obs_bool=true;
        }
        else if (line=="#SNPs") {
            std::getline (inputfile,line);
            m=atoi(line.c_str());
            SNP_bool=true;
        }
        else if (line=="#FixedEffects") {
            std::getline (inputfile,line);
            t=atoi(line.c_str());
            fixed_bool=true;
        }
        else if (line=="#FileFixedEffects") {
            std::getline (inputfile,line);
            line.copy(fXn,100);
            fixedfile_bool=true;
        }
        else if (line=="#PathGeno") {
            std::getline (inputfile,line);
            line.copy(Sd,100);
            genopath_bool=true;
        }
        else if (line=="#PathPheno") {
            std::getline (inputfile,line);
            line.copy(pd,100);
            phenopath_bool=true;
        }
        else if (line=="#TestPath") {
            std::getline (inputfile,line);
            line.copy(Ts,100);
            testpath_bool=true;
        }
        else if (line=="#DataFile") {
            std::getline (inputfile,line);
            line.copy(fDn,100);
            datafile_bool=true;
        }
        else if (line=="#DataFileHDF5") {
            std::getline (inputfile,line);
            dh5=atoi(line.c_str());
            h5_bool=true;
        }
        else if (line=="#KeepCopyOfCMatrix") {
            std::getline (inputfile,line);
            Ccop=atoi(line.c_str());
        }
        else if (line=="#TestFile") {
            std::getline (inputfile,line);
            line.copy(Ts,100);
            testfile_bool=true;
        }
        else if(line=="#MaximumIterations") {
            std::getline (inputfile,line);
            mi=atoi(line.c_str());
            maxit_bool=true;
        }
        else if (line=="#Lambda") {
            std::getline (inputfile,line);
            l=atof(line.c_str());
            lam_bool=true;
        }
        else if (line=="#Epsilon") {
            std::getline (inputfile,line);
            e=atof(line.c_str());
            eps_bool=true;
        }
        else if (line=="#BlockSize") {
            std::getline (inputfile,line);
            b=atoi(line.c_str());
            blocksize_bool=true;
        }
        else if (line=="#TestSamples") {
            std::getline (inputfile,line);
            nt=atoi(line.c_str());
            ntest_bool=true;
        }
        else if (line[0]=='/' || line.size()==0) {}
        else {
            printf("Unknown parameter in inputfile, the following line was ignored: \n");
            printf("%s\n",line.c_str());
        }
    }
    if(obs_bool) {
        if(SNP_bool) {
            if(fixed_bool) {
                if(datafile_bool) {
                    if(fixedfile_bool) {
                        if(h5_bool) {
                            if(*pst==0 && *(pst+1)==0) {
                                printf("number of observations:   \t %d\n", n);
                                printf("number of SNP effects:    \t %d\n", m);
                                printf("number of fixed effects:  \t %d\n", t);
                                printf("filename of dataset:      \t %s\n", fDn);
                                printf("filename of fixed effects:\t %s\n", fXn);
                            }
                        }
                        else {
                            printf("ERROR: filetype of dataset was not in input file or not read correctly\n");
                            return -1;
                        }
                    }
                    else {
                        printf("ERROR: filename of fixed effects was not in input file or not read correctly\n");
                        return -1;
                    }
                }
                else {
                    printf("ERROR: filename of SNP effects was not in input file or not read correctly\n");
                    return -1;
                }
            }
            else {
                printf("ERROR: number of fixed effects was not in input file or not read correctly\n");
                return -1;
            }
        }
        else {
            printf("ERROR: number of SNP effects was not in input file or not read correctly\n");
            return -1;
        }
    }
    else {
        printf("ERROR: number of observations was not in input file or not read correctly\n");
        return -1;
    }
    if(*pst==0 && *(pst+1)==0) {
        if (dh5) {
            if (genopath_bool) {
                if(phenopath_bool) {
		    printf("Dataset file is an HDF5-file \n");
                    printf("path for genotypes in dataset:      \t %s\n", Sd);
                    printf("path for phenotypes in dataset:     \t %s\n", pd);
                }
                else {
                    printf("ERROR: path for phenotypes of dataset was not in input file or not read correctly\n");
                    return -1;
                }
            }
            else {
                printf("ERROR: path for genotypes of dataset was not in input file or not read correctly\n");
                return -1;
            }
        }
        else
	  printf("Dataset file is a binary file, with as first column the phenotypical score\n");
	if(Ccop)
	  printf("A copy of the coefficient matrix will be stored throughout the computations\n");
	else
	  printf("The coefficient matrix will be read in at the beginning of every iteration to save memory\n");
        if(blocksize_bool) {
            printf("Blocksize of %d was used to distribute matrices across processes\n", b);
        }
        else {
            printf("Default blocksize of %d was used to distribute matrices across processes\n", b);
        }
        if(lam_bool)
            printf("Start value of %g was used to estimate variance component lambda\n", l);
        else
            printf("Default start value of %g was used to estimate variance component lambda\n", l);
        if(eps_bool)
            printf("Convergence criterium of %g was used to estimate variance component lambda\n", e);
        else
            printf("Default convergence criterium of %g was used to estimate variance component lambda\n", e);
        if(maxit_bool)
            printf("Maximum number of REML iterations : %d\n", mi);
        else
            printf("Default maximum number of REML iterations : %d\n", mi);
        if(dh5==0 && testfile_bool) {
            printf("Cross-validation will be performed with test set in file: \t%s\n", Ts);
            if(ntest_bool)
                printf("Cross-validation will be performed on sample with size: \t%d\n", nt);
            else {
                printf("ERROR: Number of test samples is required when cross-validation is performed\n");
                return -1;
            }
        }
        else if (dh5 && testpath_bool) {
            printf("Cross-validation will be performed with test set in path: \t%s\n", Ts);
            if(ntest_bool)
                printf("Cross-validation will be performed on sample with size: \t%d\n", nt);
            else {
                printf("ERROR: Number of test samples is required when cross-validation is performed\n");
                return -1;
            }
        }

        else
            printf("No cross-validation is performed\n");
    }
    else {
        if(testfile_bool && !ntest_bool)
            return -1;
        if(dh5 && (!genopath_bool || !phenopath_bool))
            return -1;
    }
    return 0;
}
