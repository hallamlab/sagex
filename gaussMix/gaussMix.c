
#include "gaussMix.h" 

// calculates the density of a multivariate gaussian 
// x : the point at which density is valued 
// d : the dimension of the gaussian 
// mu : mean vector 
// sig : a PSD col-major order d X d matrix 
// pi : a precomputed value for pi 
// out : space for an individual double for output 
// eps : convergence term 
// threads : number of available POSIX threads 
void mvdnorm ( double *x , int *d , double *mu , double* sig , double *pi , double *out , double *eps , int *threads ) 
{
	double dt ; 
	det( sig , d , &dt ) ; 
// fprintf( stderr , "DEBUG %e\n" , dt ) ;  
	double *prec = (double*) malloc( *d * *d * sizeof(double) ) ; 
	invPsd( sig , d , prec , eps , threads ) ;  
	int i , j ; 
	double quadForm = 0.0 ; 
	for( i = 0 ; i < *d ; i++ ) 
	{
		for( j = 0 ; j < *d ; j++ ) 
			quadForm += ( x[i] - mu[i] ) * prec[ i + *d * j ] * ( x[j] - mu[j] ) ; 
	}
	*out = -0.5 * ((double) *d) * log(2.0 * *pi) ; 
	*out += -0.5 * log(dt) - 0.5 * quadForm ;  
//	*out = exp( *out ) ; // always calculate on log scale  
	free(prec) ; 
}

// solves for a matrix of j-th dist given i-th datum 
// x : d X n col-major order data matrix 
// p : k-length probability vector 
// mu : d X k col-major matrix of means, j-th model per column 
// sig : d X d X k col-major order matrix of PSD matrices 
// pi : a precomputed value for pi 
// out : space for k X n doubles 
// eps : convergence term 
// threads : number of available POSIX threads 
void pGivenX ( double *x , int *n , int *d , int *k , double *p , double *mu , double *sig , double *pi , double *out , double *eps , int *threads ) 
{
	int i , j ; 
	double tmp , total ; 
	for( i = 0 ; i < *n ; i++ ) 
	{
		total = 0.0 ; 
		for( j = 0 ; j < *k ; j++ ) 
		{
			mvdnorm( &x[ *d * i ] , d , &mu[ *d * j ] , &sig[ *d * *d * j ] , pi , &tmp , eps , threads ) ; 
			out[ j + *k * i ] = log(p[j]) + tmp ; 
			if( j == 0 ) // log(a + b) = log(a) + log( 1 + exp( log(b) - log(a) ) ) , values are calculated in log scale to avoid over/underflows 
				total = out[ j + *k * i ] ; 
			else
				total = total + log( 1.0 + exp( out[ j + *k * i ] - total ) ) ; 
		}
		for( j = 0 ; j < *k ; j++ ) 
			out[ j + *k * i ] = exp( out[ j + *k * i ] - total ) ; 
	}
}

// Expected log likelihood 
// x : d X n col-major data matrix 
// p : k-length probability array 
// mu : d-length probability array 
// sig : d X d col-major order PSD matrix 
// pi : a precomputed value for pi 
// out : space for a single double as output  
// eps : convergence term 
// threads : number of available POSIX threads 
void eLogLik ( double *x , int *n , int *d , int *k , double *p , double *mu , double *sig , double *pmat , double *pi , double *out , double *eps , int *threads ) 
{
	//double *pmat = (double*) malloc( *k * *n * sizeof(double) ) ; 
	//pGivenX ( x , n , d , k , p , mu , sig , pi , pmat , eps , threads ) ; 
	
	int i , j ; 
	*out = 0.0 ; 
	double tmp ; 
	for( i = 0 ; i < *n ; i++ ) 
	{
		for( j = 0 ; j < *k ; j++ ) 
		{
			mvdnorm( &x[ *d * i ] , d , &mu[ *d * j ] , &sig[ *d * *d * j ] , pi , &tmp , eps , threads ) ; 
			*out += pmat[ j + *k * i ] * tmp ; 
		}
	} 
	
	//free( pmat ) ; 
}

// calculates an empirical covariance matrix 
// mat : d X n col-major data matrix 
void covMatrix ( double *mat , int *n , int *dims , double *mean , double *cov ) 
{
	int i , j , k ; 
	for( i = 0 ; i < *dims ; i++ ) // Calculate mean 
	{
		mean[i] = 0.0 ; 
		for( j = 0 ; j < *n ; j++ ) 
			mean[i] += mat[ i + *dims * j ] ; 
		mean[i] = mean[i] / ((double) *n ) ; 
	} 
	for( i = 0 ; i < *dims ; i++ ) // calc covariance matrix 
	{
		for( j = 0 ; j < *dims ; j++ ) 
		{
			cov[ i + *dims * j ] = 0.0 ; 
			for( k = 0 ; k < *n ; k++ ) 
				cov[ i + *dims * j ] += ( mat[i + *dims * k] - mean[i] ) * ( mat[j + *dims * k] - mean[j] ) ; 
			cov[ i + *dims * j ] = cov[ i + *dims * j ] / ((double) *n - 1 ) ;  
		}
	}
}

// generates unit normals filling an n X m col-major matrix 
void boxMuller( int *n , int *m , double *pi , double *out ) 
{
	int i , t ; 
	double r1 , r2 ; 
	if( *n * *m % 2 != 0 ) 
		t = *n * *m - 1 ; 
	else
		t = *n * *m ; 
	for( i = 0 ; 2*i < t ; i++ )
	{
		r1 = ((double) rand() ) / ((double) RAND_MAX) ; 
		r2 = ((double) rand() ) / ((double) RAND_MAX) ; 
		out[ 2*i ] = sqrt( -2.0 * log(r1) ) * cos( 2.0 * *pi * r2 ) ; 
		out[ 2*i + 1 ] = sqrt( -2.0 * log(r1) ) * sin( 2.0 * *pi * r2 ) ; 
	} 
	if( *n * *m % 2 != 0 ) 
	{
		r1 = ((double) rand() ) / ((double) RAND_MAX) ; 
		r2 = ((double) rand() ) / ((double) RAND_MAX) ; 
		out[ *n * *m - 1 ] = sqrt( -2.0 * log(r1) ) * cos( 20 * *pi * r2 ) ; 
	}
}

// chooses initial values for parameters 
void init ( double *x , int *n , int *d , int *k , double *p , double *eps , int* threads , double *pi ,  double *mu , double *sig , int* maxIter ) 
{
	double *sig0 = (double*) malloc( *d * *d * sizeof(double) ) ; 
	double *tmpMu = (double*) malloc( *d * sizeof(double) ) ; 
	covMatrix( x , n , d , tmpMu , sig0 ) ; 
	
	double *eigVals = (double*) malloc( *d * sizeof(double) ) ; 
	double *eigVecs = (double*) malloc( *d * *d * sizeof(double) ) ; 
	// powerIteration( sig0 , d , d , eps , eigVals , eigVecs , threads ) ; 
	double eps2 = *eps * 0.1 ; 
	psdEig ( sig0 , d , &eps2 , eigVecs , eigVals ) ; // TODO remove threads arg if no longer used 
	double *tmpMat1 = (double*) malloc( *d * *d * sizeof(double) ) ; 
	int i , j ; 
	for( i = 0 ; i < *d ; i++ ) 
	{
		for( j = 0 ; j < *d ; j++ ) 
		{
			if( i == j ) 
				tmpMat1[ i + *d * j ] = eigVals[i] / 2.0 ; 
			else
				tmpMat1[ i + *d * j ] = 0.0 ; 
		}
	}
	
	double *tmpMat2 = (double*) malloc( *d * *d * sizeof(double) ) ; 
	matTrProd ( eigVecs , tmpMat1 , d , d , d , tmpMat2 ) ; 
	matProd( tmpMat2 , eigVecs , d , d , d , tmpMat1 ) ; 
	for( i = 0 ; i < *k ; i++ ) 
		memcpy( &sig[ *d * *d * i ] , sig0 , *d * *d * sizeof(double) ) ; 
	
	// choose random initializations 
	for( i = 0 ; i < *d ; i++ ) 
	{
		for( j = 0 ; j < *d ; j++ ) 
		{
			if( i == j ) 
				tmpMat1[ i + *d * j ] = sqrt( eigVals[i] ) ; 
			else
				tmpMat1[ i + *d * j ] = 0.0 ; 
		}
	}
	matTrProd( eigVecs , tmpMat1 , d , d , d , tmpMat2 ) ; 
	matProd( tmpMat2 , eigVecs , d , d , d , tmpMat1 ) ; // tmpMat1 is now \Sigma^{1/2} for simulation 
	
	boxMuller( d , k , pi , mu ) ; // fill mu with unit normals 
	
	int one = 1 ; 
	for( i = 0 ; i < *k ; i++ ) // transform normal deviates 
	{
		matProd( tmpMat1 , &mu[ *d * i ] , d , d , &one , tmpMat2 ) ; 
		for( j = 0 ; j < *d ; j++ ) 
			mu[ j + *d * i ] = tmpMat2[ j ] + tmpMu[ j ] ; 
	}
	
	// finally, divide the probability vector into equal parts 
	for( i = 0 ; i < *k ; i++ ) 
		p[i] = 1.0 / ((double) *k) ; 
	
	// run kmeans 
	kmeans ( x , d , n , k , p , mu , sig , eps , maxIter ) ; 
	
	free( tmpMu ) ; 
	free( tmpMat1 ) ; 
	free( tmpMat2 ) ; 
	free( eigVals ) ; 
	free( eigVecs ) ; 
	free( sig0 ) ; 
}

void nextP ( double *pMat , int *n , int *d , int *k , double *p ) 
{
	int i , j ; 
	for( i = 0 ; i < *k ; i++ ) 
	{
		p[i] = 0.0 ; 
		for( j = 0 ; j < *n ; j++ ) 
			p[i] += pMat[ i + *k * j ] ; 
		p[i] = p[i] / *n ; 
	}
}

void nextMu ( double *x , int *n , int *d , int *k , double *pMat , double *out ) 
{
	double total ; 
	int i , j , l ; 
	for( j = 0 ; j < *k ; j++ ) 
	{
		total = 0.0 ; 
		for( i = 0 ; i < *n ; i++ ) 
			total += pMat[ j + *k * i ] ; 
		
		for( l = 0 ; l < *d ; l++ ) 
		{
			out[ l + *d * j ] = 0.0 ; 
			for( i = 0 ; i < *n ; i++ ) 
				out[ l + *d * j ] += x[ l + *d * i ] * pMat[ j + *k * i ] ; 
			out[ l + *d * j ] = out[ l + *d * j ] / total ; 
		}
		
	}	
}

/*
// returns a * t(b) 
// a : col-major p X q matrix 
// b : col-major r X q matrix 
// out : space for p X r doubles 
void matProdTr ( double *a , double *b , int *p , int *q , int *r , double *out ) 
{
	int i , j , k ; 
	for( i = 0 ; i < *p ; i++ ) 
	{
		for( j = 0 ; j < *r ; j++ ) 
		{
			out[ i + *p * j ] = 0.0 ; 
			for( k = 0 ; k < *q ; k++ ) 
				out[ i + *p * j ] += a[ i + *p * k ] * b[ j + *r * k ] ; 
		}
	} 
}
*/

void nextSig ( double *x , int *n , int *d , int *k , double *pMat , double *mu , double *out ) 
{
	int i , j , l , h ; 
	int one = 1 ; 
	double denom ; 
	double *tmpMat = (double*) malloc( *d * *d * sizeof(double) ) ; 
	for( j = 0 ; j < *k ; j++ ) // cycle through sigmas 
	{
		// calculate estimate denominator 
		denom = 0.0 ; 
		for( i = 0 ; i < *n ; i++ ) 
			denom += pMat[ j + *k * i ] ; 
		
		// set matrix estimate to zero 
		for( l = 0 ; l < *d ; l++ ) 
		{
			for( h = 0 ; h < *d ; h++ ) 
				out[ l + *d * h + *d * *d * j ] = 0.0 ; 
		}
		
		// calculate the estimate  
		for( i = 0 ; i < (*n) ; i++ ) // sum over samples 
		{
			for( l = 0 ; l < *d ; l++ ) 
			{
				for( h = 0 ; h < *d ; h++ ) 
					out[ l + *d * h + *d * *d * j ] += pMat[ j + *k * i ] * ( x[ l + *d * i] - mu[ l + *d * j ] ) * ( x[ h + *d * i] - mu[ h + *d * j ] ) / denom ; 
			}
		}
	} 
	free( tmpMat ) ; 
}

// attempts to initialize a gaussian mixture model estimation proceedure for a fixed value of k 
// returns -1 if failed, 0 if successful  
int gmmInit ( double *x , int *n , int *d , int *k , double *p , double *eps , int* threads ,  double *mu , double *sig , int* maxIter ) 
{ 
	double pi = atan(1.0) * 4.0 ; 
	int fails = 1 ; 
	int i , j ; 
	double determinant , smallestDeterminant ; 
	for( i = 0 ; i < *maxIter && fails > 0 ; i++ ) // reattempt initialization  
	{ 
		init ( x , n , d , k , p , eps , threads , &pi , mu , sig , maxIter ) ; 
		fails = 0 ; 
		smallestDeterminant = 1.0 ; 
		for( j = 0 ; j < *k ; j++ ) 
		{ 
			if( p[j] < *eps ) 
				fails++ ; // kmeans failure 
			det( &sig[ *d * *d * i ] , d , &determinant ) ; 
			if( determinant < smallestDeterminant ) 
				smallestDeterminant = determinant ; 
		} 
		if( smallestDeterminant < *eps ) 
			fails++ ; // kmeans failure 
	} 
	if( fails > 0 ) 
		return -1 ; 
	else 
		return 0 ; 
} 

// returns -1 if failed, 0 if successful 
int fitMixture( double *x , int *n , int *d , int *k , double *eps , double *p , double *mu , double *sig , int *maxIter , int* threads ) 
{
	double pi = atan(1.0) * 4.0 ; 
	
	/*
// fprintf( stderr , "DEBUG: initializing\n" ) ; 
	// initialize parameters 
	int fails = 1 ; 
	int i , j ;  
	double determinant , smallestDeterminant ; 
	for( i = 0 ; i < *maxIter && fails > 0 ; i++ ) // retry initialization  
	{
		init ( x , n , d , k , p , eps , threads , &pi , mu , sig , maxIter ) ; 
		zeros = 0 ; 
		smallestDeterminant = 1 ; 
		for( j = 0 ; j < *k ; j++ ) 
		{
			if( p[j] < *eps ) 
				fails++ ; // kmeans failure 
			det( &sig[ *d * *d * i ] , d , &determinant ) ; 
			if( determinant < smallestDeterminant ) 
				smallestDeterminant = determinant ; 
		}
		if( smallestDeterminant < *eps ) 
			fails++ ; // kmeans failure 
	}
	if( fails > 0 ) 
	{
		fprintf( stderr , "ERROR: kmeans failure\n" ) ; 
		return ; 
	} 
// fprintf( stderr , "DEBUG: first lik\n" ) ;
	
	double determinant ; 
	double smallestDeterminant = 1 ; 
	for( i = 0 ; i < *k ; i++ ) 
	{ 
		det( &sig[ *d * *d * i ] , d , &determinant ) ; 
		if( determinant < smallestDeterminant ) 
			smallestDeterminant = determinant ; 
	}
	if( smallestDeterminant < *eps ) 
		fprintf( stderr , "WARNING: Covariance numerically singular (determinant = %e). Reduce argument -k\n" , smallestDeterminant ) ; 
	*/
	
	// calculate first value for log likelihood 
	double prevLik , lik ; 
// fprintf( stderr , "calculating first eloglik\n" ) ; 
	// eLogLik ( x , n , d , k , p , mu , sig , &pi , &prevLik , eps , threads ) ; 
	
	// initialize variables for main loop 
	double *tmpP = (double*) malloc( *k * sizeof(double) ) ; 
	double *pMat = (double*) malloc( *k * *n * sizeof(double) ) ; 
	double *tmpMu = (double*) malloc( (*d) * (*k) * sizeof(double) ) ; 
	double *tmpSig = (double*) malloc( *d * *d * *k * sizeof(double) ) ; 
	
	pGivenX ( x , n , d , k , p , mu , sig , &pi , pMat , eps , threads ) ; 
	nextP ( pMat , n , d , k , tmpP ) ; 
	nextMu ( x , n , d , k , pMat , tmpMu ) ; 
	nextSig ( x , n , d , k , pMat , tmpMu , tmpSig ) ; 
// fprintf( stderr , "DEBUG tmpMu: " ) ; for( i = 0 ; i < *d * *k ; i++ ) fprintf( stderr , "%e " , tmpMu[i] ) ; fprintf( stderr , "\n" ) ; 
	
	eLogLik ( x , n , d , k , tmpP , tmpMu , tmpSig , pMat , &pi , &prevLik , eps , threads ) ; 
	
	memcpy( p , tmpP , *k * sizeof(double) ) ; 
	memcpy( mu , tmpMu , (*d) * (*k) * sizeof(double) ) ;  
	memcpy( sig , tmpSig , *d * *d * *k * sizeof(double) ) ; 
	
// fprintf( stderr , "DEBUG: main loop\n" ) ;
	double err = *eps + 1.0 ; 
	int contin = 1 ; 
	int iter = 0 ;  
// fprintf( stderr , "running main loop\n" ) ; 
	while( (err > *eps && contin > 0) && iter < *maxIter && contin > 0 ) 
	{
// fprintf( stderr , "DEBUG: err: %e, eps: %e, iter: %i\n" , err , *eps , iter ) ;
		iter += 1 ; 
		
		pGivenX ( x , n , d , k , p , mu , sig , &pi , pMat , eps , threads ) ; 
		nextP ( pMat , n , d , k , tmpP ) ; 
		nextMu ( x , n , d , k , pMat , tmpMu ) ; 
		nextSig ( x , n , d , k , pMat , tmpMu , tmpSig ) ; 
// fprintf( stderr , "DEBUG tmpMu: " ) ; for( i = 0 ; i < *d * *k ; i++ ) fprintf( stderr , "%e " , tmpMu[i] ) ; fprintf( stderr , "\n" ) ;
		
		eLogLik ( x , n , d , k , tmpP , tmpMu , tmpSig , pMat , &pi , &lik , eps , threads ) ; 
		if( lik == lik ) // Check for nans 
		{
			if( 1 > 0 ) // ( lik > prevLik ) TODO   
			{ 
// fprintf( stderr , "DEBUG: lik - prevLik: %e\n" , lik - prevLik ) ; 
				// prevLik = lik ; 
				memcpy( p , tmpP , *k * sizeof(double) ) ; 
				memcpy( mu , tmpMu , (*d) * (*k) * sizeof(double) ) ; 
				memcpy( sig , tmpSig , *d * *d * *k * sizeof(double) ) ; 
			} 
			else
			{ 
				fprintf( stderr , "DEBUG: lik !> prevLik, %e !> %e on iter %i\n" , lik , prevLik , iter ) ; 
				contin = -2 ; 	
			} 
			err = fabs( lik - prevLik ) ; 
			prevLik = lik ; 
			if( err < *eps ) 
				contin = 0 ; 
		}
		else
		{ 
			fprintf( stderr , "DEBUG: lik == nan\n" ) ; 
			contin = -1 ; 
		} 
	}
// fprintf( stderr , "contin: %i\n" , contin ) ; 
	
	// free memory 
	free( tmpP ) ; 
	free( tmpMu ) ; 
	free( tmpSig ) ; 
	free( pMat ) ; 
	
	if( contin < 0 ) 
		return -1 ; 
	else 
		return 0 ; 

// fprintf( stderr , "DEBUG:\n" ) ; 
// fprintf( stderr , "p: " ) ; for( i = 0 ; i < *k ; i++ ) fprintf( stderr , "%e " , p[i] ) ; fprintf( stderr , "\n" ) ; 
// for( i = 0 ; i < *k ; i++ ) { fprintf( stderr , "mu %i : " , i ) ; for( j = 0 ; j < *d ; j++ ) fprintf( stderr , "%e " , mu[j + *d * i] ) ; fprintf( stderr , "\n" ) ; }  int kk ; 
// for( kk = 0 ; kk < *k ; kk++ ){ fprintf( stderr , "sig %i:\n" , kk ) ; for( i = 0 ; i < *d ; i++ ){ for( j = 0 ; j < *d ; j++ ) fprintf( stderr , "%e " , sig[ i + *d * j + *d * *d * kk ] ) ; fprintf( stderr , "\n" ) ; } }
}













