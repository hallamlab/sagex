#ifndef HELPER_H
#define HELPER_H

#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <getopt.h>

struct Options {
    /* Input files for SAGEX */
    char *input; //The fasta formatted file containing the SAG
    char *Gm; //The fasta formatted metagenome
    char *output; //The output file for the contigs deemed a part of the SAG's k-mer space
    char *kmerPCA; //Optional output file for k-mer PCA
    char *kmerFreq; //Optional output file for the SAG tetranucleotide frequency vectors

    /* Input settings*/
    int k; // The number of Gaussians to use
    int chopSize;
    int overlap;
    double Alpha;
    double Beta;
    double eps;
    double mahalaMultiple; //The Mahalanobis multiple
    int itLen;
    int itMax;
    int lcsCut;
    int threads; //Number of threads to use
    int seed; //The seeding number for random number generation
    int errFlag;
    bool proportionFlag; //k-mer frequencies will be based on tetranucleotide proportions of whole contigs
    bool hFlag; //Flag triggering the Print_help function
    bool fixK; //Flag indicating to roll back the number of Gaussians if desired number cannot be used
    bool verbose; //Self explanatory
    bool problem; //For getopt

    /*  Constructor */
    Options() {
        input = (char *) malloc (100 * sizeof(char)); //TODO:change this to a more informative variable
        Gm = (char *) malloc (100 * sizeof(char));
        output = (char *) malloc (100 * sizeof(char));
        output[0] = '\0';
        kmerFreq = (char *) malloc (100 * sizeof(char));
        kmerFreq[0] = '\0';
        kmerPCA = (char *) malloc (100 * sizeof(char));
        kmerPCA[0] = '\0';

        chopSize = 0;
        overlap = 0;
        k = 2;
        Alpha = 0.05;
        Beta = 0.8;
        eps = 0.0000001;
        mahalaMultiple = 1.0 ;
        itLen = 2000;
        itMax = 10000;
        lcsCut = -1 ;
        threads = 1;
        seed = -1 ;

        proportionFlag = true;
        hFlag = false;
        fixK = false;
        verbose = false;
        problem = false;

        errFlag = opterr = 0;
    };

    /* Destructor */
    ~Options() {
        ;
    }

    /* Functions */
    void Print_usage(char* arg);
    void Print_help();
    void ReviewArguments(char* argv[]);
    int FetchArguments(int argc, char* argv[]);
};

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);

std::vector<std::string> split(const std::string &s, char delim);

#endif
