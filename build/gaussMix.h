
#include <math.h> 
#include "matrix.h" 
#include "kmeans.h"

// fits a gaussian mixture model 
// x : d X n col-major data matrix 
// n : number of data 
// d : dimension of data 
// k : desired number of gaussians to fit in mixture 
// eps : eigenvector calculation convergence term 
// p : output, k-length probability vector 
// mu : output, d X k mean vetors 
// sig : output, d X d X k covaraince matrices 
// maxIter : maximum number of iterations allowed in estimation 
// threads : number of POSIX threads available 
void fitMixture( double *x , int *n , int *d , int *k , double *eps , double *p , double *mu , double *sig , int *maxIter , int* threads ) ; 

