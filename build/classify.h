
#include "matrix.h" 
#include "gammaDist.h" 
// #include "shared.h" 
#include "chopCounter.h" 
#include "gaussMix.h"
#include "gmmPval.h"

// Classifies contigs as IN a SAG 
// sag : list of standard-length sag sequences 
// sagN : length of the sag array 
// names : integer titles equating only when a metagenomic sequence is from the same contig 
// gm : list of standard-length sag / metbag sequences 
// alpha : type 1 error rate for classifying an element of gm as IN the SAG 
// beta : proportion of shared names that need to be IN the SAG 
// threads : the desired number of posix threads for this task 
// eps : PCA convergence error parameter 
// minIter : PCA minimum number of iterations 
// maxIter : PCA maximum number of iterations 
// chopsSize : 
// overlap : 
// k : number of gaussians in the gaussian mixture model modelling the SAG 
// verbose : describe process in stderr if > 0 
// out : a pointer to be allocated with the int-names of contigs which have been classified as IN the SAG 
void classify ( char **sag , int sagN , char **gm , int gmN , double alpha , double beta , int threads , double eps , int minIter , int maxIter , int chopSize , int overlap , int k , int verbose , int **out ) ; 

// Standardizes the columns of a matrix 
void colStandardize( double *mat , int *rows , int *cols , double *out ) ; 

// out : enough space for dims X dims doubles 
void corrMat ( double *mat , int *n , int *dims , double *out ) ; 

void covMat ( double *mat , int *n , int *dims , double *mean , double *cov ) ; 

