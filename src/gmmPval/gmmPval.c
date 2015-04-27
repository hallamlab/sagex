
#include "gmmPval.h"

// mu : dim - vector  
// invSqrtSig : dim X dim vector 
double calcSingleStat ( double *x , int *dim , double *mu , double *invSqrtSig ) 
{ 
	double *tmpVec1 = (double*) malloc( *dim * sizeof(double) ) ; // TODO move malloc out of loops 
	double *tmpVec2 = (double*) malloc( *dim * sizeof(double) ) ; 
	int i ; 
	for( i = 0 ; i < *dim ; i++ ) 
		tmpVec1[i] = x[i] - mu[i] ; 
	
	int j ; 
	for( i = 0 ; i < *dim ; i++ ) 
	{ 
		tmpVec2[i] = 0.0 ; 
		for( j = 0 ; j < *dim ; j++ ) 
			tmpVec2[i] += tmpVec1[j] * invSqrtSig[ j + *dim * i ] ; 
	} 
	
	double out = 0.0 ; 
	for( i = 0 ; i < *dim ; i++ ) 
		out += tmpVec2[i] * tmpVec2[i] ; 
	
	free( tmpVec1 ) ; 
	free( tmpVec2 ) ; 
	
	return out ; 
} 

// Calculates the minimum Malahanobis distance from x to a centroid  
// k : number of gaussians 
// mu : dim X k 
// invSqrtSig : dim X dim X k matrix of inverse square-root matrices  
double calcStat ( double *x , int *k , int *dim , double *mu , double *invSqrtSig ) 
{ 
	double out = calcSingleStat ( x , dim , mu , invSqrtSig ) ; 
	double tmp ; 
	int i ; 
	for( i = 1 ; i < *k ; i++ ) 
	{ 
		tmp = calcSingleStat ( x , dim , &mu[ *dim * i ] , &invSqrtSig[ *dim * *dim * i ] ) ; 
		if( tmp < out ) 
			out = tmp ; 
	} 
	return out ; 
} 

// generate univariate normals 
// pi : 3.14159... etc 
// out : space for two doubles  
void boxMuller ( double *pi , double *out ) 
{
	double u1 = ((double) rand() ) / ((double) RAND_MAX ) ; 
	double u2 = ((double) rand() ) / ((double) RAND_MAX ) ; 
	
	out[0] = sqrt( -2.0 * log(u1) ) * cos( 2.0 * *pi * u2 ) ; 
	out[1] = sqrt( -2.0 * log(u1) ) * sin( 2.0 * *pi * u2 ) ; 
}

// fills a matrix with univariate normals 
// out : an n X m matrix 
void matMuller ( int *n , int *m , double *out ) 
{
	double pi = 4.0 * atan(1.0) ; 
	int t = *n * *m ; 
	int i ; 
	for( i = 0 ; 2*i+1 < t ; i++ ) 
		boxMuller ( &pi , &out[2*i] ) ; 
	if( t % 2 != 0 && t > 0 ) 
	{
		double tmp[2] ; 
		boxMuller ( &pi , &tmp[0] ) ; 
		out[ t-1 ] = tmp[0] ; 
	}
}

// simulates values from a GMM 
// n : the number of data to simulate  
// d : the dimension 
// k : the number of gaussians 
// p : probability vector for choice variable 
// mu : d X k col-major order matrix of means 
// sig : d X d X k col-major order matrix of covariances 
// eps : tolerance for eigen-vector calculations 
// threads : number of POSIX threads available 
// out : a space for d X n doubles 
void simulateGMM ( int *n , int *d , int *k , double *p , double *mu , double *sig , double *eps , int *threads , double *out ) 
{
	double *sqrtSig = (double*) malloc( *d * *d * *k * sizeof(double) ) ; 
	double *eigVecs = (double*) malloc( *d * *d * sizeof(double) ) ; 
	double *eigVals = (double*) malloc( *d * sizeof(double) ) ; 
	double *diagMat = (double*) malloc( *d * *d * sizeof(double) ) ; 
	int i , j ; 
	for( i = 0 ; i < *d ; i++ ) // set diagMat entries to zero 
	{
		for( j = 0 ; j < *d ; j++ ) 
			diagMat[ i + *d * j ] = 0.0 ; 
	}
	for( i = 0 ; i < *k ; i++ ) // calc sqrt-mats 
	{
		// powerIteration ( &sig[ *d * *d * i ] , d , d , eps , eigVals , eigVecs , threads ) ; 
		double eps2 = *eps * 0.1 ; 
		psdEig ( &sig[ *d * *d * i ] , d , &eps2 , eigVecs , eigVals ) ; // TODO remove threads arg if not used 
		for( j = 0 ; j < *d ; j++ ) 
			diagMat[ j + *d * j ] = sqrt(eigVals[j]) ; 
		matProdP ( eigVecs , diagMat , d , d , d , &sqrtSig[ *d * *d * i ] , threads ) ; 
	}
	
	// generate a distribution function for the choice variable 
	double *pp = (double*) malloc( *k * sizeof(double) ) ; 
	pp[0] = p[0] ; 
	for( i = 1 ; i < *k ; i++ ) 
		pp[i] = pp[i-1] + p[i] ;  
	
	// generate unit normal simulants 
	matMuller ( d , n , out ) ; 
	
	// run them through the GMM 
	double r ; 
	int choice ; 
	double *tmpVec = (double*) malloc( *d * sizeof(double) ) ; 
	int one = 1 ; 
	for( i = 0 ; i < *n ; i++ ) 
	{
		// generate a choice 
		r = ((double) rand() ) / ((double) RAND_MAX ) ; 
		choice = -1 ; 
		if( r <= pp[0] ) 
			choice = 0 ; 
		for( j = 1 ; j < *k && choice < 0 ; j++ ) 
		{
			if( r > pp[j-1] && r <= pp[j] ) 
				choice = j ; 
		}
	
		// transform the simulant 
		matProd ( &sqrtSig[ *d * *d * choice ] , &out[ *d * i ] , d , d , &one , tmpVec ) ; 
		for( j = 0 ; j < *d ; j++ ) 
			out[j + *d * i] = tmpVec[j] + mu[ j + *d * choice ] ; 
	}
	
	free( tmpVec ) ; 
	free( sqrtSig ) ; 
	free( eigVecs ) ; 
	free( eigVals ) ; 
	free( diagMat ) ; 
	free( pp ) ; 
}

/*
void quickSort ( double *x , int *n ) 
{
	if( *n <= 1 ) 
		return ; 
	double tmp ; 
	if( *n == 2 ) 
	{
		if( x[0] > x[1] ) 
		{
			tmp = x[1] ; 
			x[1] = x[0] ; 
			x[0] = tmp ; 
		}
		return ; 
	}
	
	int pivot = *n - 1 ;  
	
	int i ; 
	int j = 0 ; 
	for( i = 0 ; i < *n-1 ; i++ ) 
	{
		if( x[i] < x[pivot] ) 
		{
			tmp = x[i] ; 
			x[i] = x[j] ; 
			x[j] = tmp ; 
			j++ ; 
		}
	}
	tmp = x[pivot] ; 
	x[pivot] = x[j] ; 
	x[j] = tmp ; 
	pivot = j ; 
	
	i = pivot  ; 
	quickSort( x , &i ) ; 
	i = *n - pivot - 1 ; 
	quickSort( &x[pivot+1] , &i ) ; 
}
*/

// parametric bootstrap calculating quantiles of minimum distance from the means 
// q : quantile value desired 
// n : bootsrap count 
// p : k-length probability vector for GMM of k gaussians 
// mu : d X k matrix of k means for the GMM 
// sig : d X d X k matrix of covaraince matricies for the GMM 
// eps : eigen-vector calculation tolerance 
// threads : nubmer of POSIX threads available for calculation 
// out : space for a single double 
void getGMMQuantile ( double *q , double *p , double *mu , double *sig , int *n , int *d , int *k , double *eps , int *threads, double *out ) 
{
	// simulate bootstraps from GMM 
	double *sims = (double*) malloc( *d * *n * sizeof(double) ) ; 
	double *stats = (double*) malloc( *n * sizeof(double) ) ; 
	simulateGMM ( n , d , k , p , mu , sig , eps , threads , sims ) ; 
	
	// calculate inverse square root matrices 
	int i , j , l ; 
	double *invSqrtSig = (double*) malloc( *d * *d * *k * sizeof(double) ) ; // todo : free this  
	double *tmpMat = (double*) malloc( *d * *d * sizeof(double) ) ; // todo : free this 
	double *tmpVec = (double*) malloc( *d * sizeof(double) ) ; // todo : free this 
	for( i = 0 ; i < *k ; i++ ) 
	{ 
		psdEig ( &sig[ *d * *d * i ] , d , eps , tmpMat , tmpVec ) ; 
		for( j = 0 ; j < *d ; j++ ) 
		{ 
			for( l = 0 ; l < *d ; l++ ) 
				invSqrtSig[ j + *d * l + *d * *d * i ] = tmpMat[ j + *d * l ] / sqrt(tmpVec[l]) ; 
		} 
	} 
	
	// calculate minimum distances  
	// double dist, min ; 
	for( i = 0 ; i < *n ; i++ ) 
	{ 
		/* 
		for( j = 0 ; j < *k ; j++ ) 
		{
			dist = 0.0 ; 
			for( l = 0 ; l < *d ; l++ ) 
				dist += ( mu[l + *d * j] - sims[l + *d * i] )*( mu[l + *d * j] - sims[l + *d * i] ) ; 
			dist = sqrt(dist) ; 
			if( j == 0 ) 
				min = dist ; 
			else
			{
				if( dist < min ) 
					min = dist ; 
			}
		}
		stats[i] = min ; 
		*/ 
		stats[i] = calcStat ( &sims[ *d * i ] , k , d , mu , invSqrtSig ) ; 
	}
	
	// sort 
	quickSort( stats , n , NULL ) ; 
	
	// get quantile 
	int m = floor( ((double) *q) * ((double) *n) ) ; 
	if( m > *n )
		m = *n ; 
	*out = stats[m] ; 
	
	free( sims ) ; 
	free( stats ) ; 
	free( invSqrtSig ) ; 
	free( tmpMat ) ; 
	free( tmpVec ) ; 
}





