
#include <math.h> 
#include "matrix.h" 

// Calculates the minimum Malahanobis distance from x to a centroid  
// k : number of gaussians 
// mu : dim X k 
// invSqrtSig : dim X dim X k matrix of inverse square-root matrices  
double calcStat ( double *x , int *k , int *dim , double *mu , double *invSqrtSig ) ; 

// parametric bootstrap calculating quantiles of minimum distance from the means 
// q : quantile value desired 
// n : bootsrap count 
// p : k-length probability vector for GMM of k gaussians 
// mu : d X k matrix of k means for the GMM 
// sig : d X d X k matrix of covaraince matricies for the GMM 
// eps : eigen-vector calculation tolerance 
// threads : nubmer of POSIX threads available for calculation 
// out : space for a single double 
void getGMMQuantile ( double *q , double *p , double *mu , double *sig , int *n , int *d , int *k , double *eps , int *threads, double *out ) ; 

