
//#include "classify.h"
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>

// Standardizes the columns of a matrix 
void colStandardize( double *mat , int *rows , int *cols , double *out )
{
        int i , j ;
        double mean ;
        double var ;
        for( i = 0 ; i < *cols ; i++ )
        {
                mean = 0.0 ;
                var = 0.0 ;
                for( j = 0 ; j < *rows ; j++ )
                        mean += mat[ j + *rows * i ] ;
                mean = mean / ((double) *rows) ;
                for( j = 0 ; j < *rows ; j++ )
                        var += ( mat[j + *rows * i] - mean  ) * ( mat[j + *rows * i] - mean ) ;
                var = var / ((double) (*rows - 1) ) ;

                for( j = 0 ; j < *rows ; j++ )
                        out[ j + *rows * i ] = (mat[ j + *rows * i ] - mean ) / sqrt(var) ;
        }
}

// out : enough space for dims X dims doubles 
void corrMat ( double *mat , int *n , int *dims , double *out )
{
	int i , j , k ; 
	double *m = (double*) malloc( *dims * sizeof(double) ) ; // array of means 
	double *v = (double*) malloc( *dims * sizeof(double) ) ; // array of variances 
	for( i = 0 ; i < *dims ; i++ ) // calc means 
	{
		m[i] = 0.0 ; 
		for( j = 0 ; j < *n ; j++ ) 
			m[i] += mat[ j + *n * i ] ;   
		m[i] = m[i] / ((double) *n ) ; 
	}
	for( i = 0 ; i < *dims ; i++ ) // calc vars 
	{
		v[i] = 0.0 ; 
		for( j = 0 ; j < *n ; j++ )
			v[i] += ( mat[j + *n * i] - m[i] ) * ( mat[j + *n * i] - m[i] ) ; 
		v[i] = v[i] / ((double) *n - 1 ) ; 
	}
	for( i = 0 ; i < *dims ; i++ ) // calc corrs 
	{
		for( j = 0 ; j < *dims ; j++ ) 
		{
			out[ i + *dims * j ] = 0.0 ; 
			for( k = 0 ; k < *n ; k++ ) 
				out[ i + *dims * j ] += ( mat[k + *n * i] - m[i] ) * ( mat[k + *n * j] - m[j] ) ; 
			out[ i + *dims * j ] = out[ i + *dims * j ] / ( sqrt( v[i] * v[j] ) * ((double) *n - 1 ) ) ; 
		}
	}
	free( m ) ; 
	free( v ) ; 
}

void covMat ( double *mat , int *n , int *dims , double *mean , double *cov ) 
{
	int i , j , k ; 
	for( i = 0 ; i < *dims ; i++ ) // Calculate mean 
	{
		mean[i] = 0.0 ; 
		for( j = 0 ; j < *n ; j++ ) 
			mean[i] += mat[ j + *n * i ] ; 
		mean[i] = mean[i] / ((double) *n ) ; 
	} 
	for( i = 0 ; i < *dims ; i++ ) // calc covariance matrix 
	{
		for( j = 0 ; j < *dims ; j++ ) 
		{
			cov[ i + *dims * j ] = 0.0 ; 
			for( k = 0 ; k < *n ; k++ ) 
				cov[ i + *dims * j ] += ( mat[k + *n * i] - mean[i] ) * ( mat[k + *n * j] - mean[j] ) ; 
			cov[ i + *dims * j ] = cov[ i + *dims * j ] / ((double) *n - 1 ) ;  
		}
	}
}












