
#include "kmeans.h"

// initialization scheme from k-means++ 
void initKMeans ( double *x , int *d , int *n , int *k , double *mu ) 
{
	int *idx = (int*) malloc( *k * sizeof(int) ) ; 
	double *dList = (double*) malloc( *n * sizeof(double) ) ; 
	
	idx[0] = rand() % *n ; 
	
	double tmp , total , r ; 
	int i , j , l , h , flag ; 
	for( i = 1 ; i < *k ; i++ ) // create k centroids 
	{
		for( j = 0 ; j < *n ; j++ ) // consider all points as potential centroids 
		{
			flag = 0 ; 
			for( l = 0 ; l < i && flag == 0 ; l++ ) 
			{
				if( idx[l] == j ) 
				{
					flag = 1 ; 
					dList[j] = 0.0 ; 
				}
			}
			if( flag == 0 ) // points may only be used for at most one centroid 
			{
				// initialize dList entry to first centroid 
				dList[j] = 0.0 ; 
				for( l = 0 ; l < *d ; l++ ) 
				{
					dList[j] += pow( fabs( x[ l + *d * j ] - x[ l + *d * idx[0] ] ) , 2.0 ) ; 
				}
				dList[j] = sqrt( dList[j] ) ; 
				for( l = 1 ; l < i ; l++ ) // check for nearer centroids 
				{
					tmp = 0.0 ; 
					for( h = 0 ; h < *d ; h++ ) 
					{ 
// fprintf( stderr , "h: %i, d: %i , i: %i, j: %i, l: %i, idx[l]: %i\n" , h , *d, i, j , l , idx[l] ) ; 
						tmp += pow( fabs( x[ h + *d * j] - x[ h + *d * idx[l] ] ) , 2.0 ) ; 
					} 
					if( tmp < dList[j] ) 
						dList[j] = tmp ; 
				}
			}
		}
		// pick the next centroid 
		total = 0.0 ; 
		for( j = 0 ; j < *n ; j++ ) 
			total += dList[j] * dList[j] ; 
		r = ((double) (rand() % 1000000))/1000000.0 ; 
		idx[i] = -1 ; 
		tmp = 0.0 ; 
		for( j = 0 ; j < *n && idx[i] < 0 ; j++ ) 
		{
			tmp += dList[j] * dList[j] ; 
			if( r <= tmp / total ) 
				idx[i] = j ; 
		} 
	}
	
	free( idx ) ; 
	free( dList ) ; 
}

// kmeans by Lloyd's algorithm 
// x : a col-major order d X n data matrix 
// k : number of clusters 
// p : k-length vector, output for probability of inclusion 
// mu : col-major d X k matrix of initial estimates and output for centroids 
// cov : col-major d X d X k matrix for output of covariances 
// eps : convergence tolerance 
// maxIter : maximum number of iterations allowed 
void kmeans ( double *x , int *d , int *n , int *k , double *p , double *mu , double *cov , double *eps , int *maxIter ) 
{
	initKMeans ( x , d , n , k , mu ) ; 
	
	int *memberCount = (int*) malloc( *k * sizeof(int) ) ; 
	double *muTemp = (double*) malloc( *d * *k * sizeof(double) ) ; 
	double *errList = (double*) malloc( *k * sizeof(double) ) ; 
	
	double err = *eps + 1.0 ; 
	int i , j , l , h , g , winner ; 
	double tmp ; 
	for( i = 0 ; i < *maxIter && err > *eps ; i++ ) 
	{
		memcpy( muTemp , mu , *d * *k * sizeof(double) ) ; 
		for( j = 0 ; j < *k ; j++ ) 
		{
			for( l = 0 ; l < *d ; l++ ) 
				mu[ l + *d * j ] = 0.0 ; 
			memberCount[j] = 0 ; 
		} 
		
		for( j = 0 ; j < *n ; j++ ) 
		{
			// find distances to centroids 
			for( l = 0 ; l < *k ; l++ ) 
			{
				errList[l] = 0.0 ; 
				for( g = 0 ; g < *d ; g++ ) 
					errList[l] += fabs( muTemp[ g + *d * l ] - x[ g + *d * j ] ) ; 
			}
			
			// find nearest centroid 
			winner = 0 ; 
// printf( "DEBUG. errList[0]: %e\n" , errList[0] ) ; 
			for( l = 1 ; l < *k ; l++ ) 
			{
// printf( "DEBUG. errList[%i]: %e\n" , l , errList[l] ) ; ;  
				if( errList[l] < errList[winner] ) 
					winner = l ; 
// printf( "DEBUG. winner: %i\n" , winner ) ; 
			}
			
			// assign the datum to the mean of the winner centroid 
			memberCount[ winner ] += 1 ; 
			for( g = 0 ; g < *d ; g++ ) 
				mu[ g + *d * winner ] += x[ g + *d * j ] ; 
		}
		// calculate means 
		for( j = 0 ; j < *k ; j++ ) 
		{
			if( memberCount[j] > 0 ) 
			{
				for( g = 0 ; g < *d ; g++ ) 
					mu[ g + *d * j ] = mu[ g + *d * j ] / ((double) memberCount[j]) ; 
			}
			else
			{
				for( g = 0 ; g < *d ; g++ ) 
					mu[ g + *d * j ] = 0.0 ; // error case  
			}
		}
		
		// cacluate total distances moved since last centroids 
		err = 0.0 ; 
		for( j = 0 ; j < *k ; j++ ) 
		{
			tmp = 0.0 ; 
			for( g = 0 ; g < *d ; g++ ) 
				tmp += ( mu[ g + *d * j ] - muTemp[ g + *d * j ] ) * ( mu[ g + *d * j ] - muTemp[ g + *d * j ] ) ; 
			tmp = sqrt(tmp) ; 
			err += tmp ; 
		}
// printf( "DEBUG. iter: %i, err: %e\n" , i , err ) ; 
	}
	
	// calculate covariances 
	for( l = 0 ; l < *k ; l++ ) 
	{
		memberCount[l] = 0 ; 
		for( i = 0 ; i < *d ; i++ ) 
		{
			for( j = 0 ; j < *d ; j++ ) 
				cov[ i + *d * j + *d * *d * l ] = 0.0 ; 
		}
	}
	// find membership 
	for( i = 0 ; i < *n ; i++ ) 
	{
		for( j = 0 ; j < *k ; j++ ) 
		{
			errList[j] = 0.0 ; 
			for( l = 0 ; l < *d ; l++ ) 
				errList[j] += fabs( mu[ l + *d * j ] - x[ l + *d * i ] ) ; 
		}
		winner = 0 ; 
		for( j = 1 ; j < *k ; j++ ) 
		{
			if( errList[j] < errList[winner] ) 
				winner = j ; 
		}
		memberCount[ winner ] += 1 ; 
		for( j = 0 ; j < *d ; j++ ) 
		{
			for( l = 0 ; l < *d ; l++ ) 
				cov[ j + *d * l + *d * *d * winner ] += ( x[ j + *d * i ] - mu[ j + *d * winner ] ) * ( x[ l + *d * i ] - mu[ l + *d * winner ] ) ; 
		}
	}
	// calculate covariances 
	for( i = 0 ; i < *d ; i++ ) 
	{
		for( j = 0 ; j < *d ; j++ ) 
		{
			for( l = 0 ; l < *k ; l++ ) 
			{
				if( memberCount[l] > 1 ) 
					cov[ i + *d * j + *d * *d * l ] = cov[ i + *d * j + *d * *d * l ] / ((double) memberCount[l] ) ; 
				else 
					cov[ i + *d * j + *d * *d *l ] = -1.0 ; // error case 
			}
		}
	}
	// calculate p 
	tmp = 0.0 ; 
	for( i = 0 ; i < *k ; i++ ) 
		tmp += (double) memberCount[i] ; 
	for( i = 0 ; i < *k ; i++ ) 
		p[i] = ((double) memberCount[i]) / tmp ; 
	
	free( memberCount ) ; 
	free( muTemp ) ; 
	free( errList ) ; 
}








