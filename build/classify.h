
#include "matrix.h" 
#include "gammaDist.h" 
// #include "shared.h" 
//#include "chopCounter.h" 
#include "count.h"
#include "gaussMix.h"
#include "gmmPval.h"

// Classifies contigs as IN a SAG 
// sag : list of standard-length sag sequences 
// sagN : length of the sag array 
// sagNames : a pointer to the list of sag contig names 
// names : integer titles equating only when a metagenomic sequence is from the same contig 
// gm : list of standard-length metagenome / metbag sequences 
// gmNames : a pointer to the list of metagenome contig names 
// alpha : type 1 error rate for classifying an element of gm as IN the SAG 
// beta : proportion of shared names that need to be IN the SAG 
// threads : the desired number of posix threads for this task 
// eps : PCA convergence error parameter 
// minLength : minimum length for contigs
// maxIter : PCA maximum number of iterations , max iterations in gamma dist integration or max bootstraps in qGMM integration  
// chopsSize : 
// overlap : 
// k : number of gaussians in the gaussian mixture model modelling the SAG 
// verbose : describe process in stderr if > 0 
// kmerFreq : writes kmers to file specified 
// out : a pointer to be allocated with the int-names of contigs which have been classified as IN the SAG 
void classify ( char **sag , int sagN , char **sagNames , char **gm , int gmN , char **gmNames , double alpha , double beta , int threads , double eps , int minLength , int maxIter , int chopSize , int overlap , int proportion, int k , int verbose , char *kmerFreq , char *kmerPCA, int **out, char *output) ; 

// Standardizes the columns of a matrix 
void colStandardize( double *mat , int *rows , int *cols , double *out ) ; 

// out : enough space for dims X dims doubles 
void corrMat ( double *mat , int *n , int *dims , double *out ) ; 

void covMat ( double *mat , int *n , int *dims , double *mean , double *cov ) ; 

