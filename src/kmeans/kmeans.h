
#include <stdlib.h> 
#include <stdio.h>
#include <math.h> 
#include <string.h> 

// kmeans by Lloyd's algorithm 
// x : a col-major order d X n data matrix 
// k : number of clusters 
// p : k-length vector, output for probability of inclusion 
// mu : col-major d X k matrix of initial estimates and output for centroids 
// cov : col-major d X d X k matrix for output of covariances 
// eps : convergence tolerance 
// maxIter : maximum number of iterations allowed 
void kmeans ( double *x , int *d , int *n , int *k , double *p , double *mu , double *cov , double *eps , int *maxIter ) ; 

