
#include <stdlib.h> 
#include <stdio.h> 
#include <limits.h> 
#include <string.h> 

#include "skey.c" 

int main( int argc , char **argv ) 
{
	if( argc < 4 ) 
	{ 
		printf( "Please provide two strings (arg 1 & 2) and kmer length (arg3)\n" ) ; 
		return -1 ; 
	} 
	
	int k = atoi( argv[3] ) ; 
	if( k < 1 ) 
	{
		printf( "kmer lengths must be > 0\n" ) ; 
		return -2 ; 
	}
	
	int N = strlen(argv[1]) ; 
	
	int n , preN , lastN ; 
	sHashGetKeySize ( k , &n , &preN , &lastN ) ; 
	printf( "k: %i, n: %i, preN: %i, lastN: %i\n" , k , n , preN , lastN ) ; 
	
	ull *hash = (ull*) malloc( n * (N-k+1) * sizeof(ull) ) ; 
	printf( "keys to store: %i\n" , N-k+1 ) ; 
	
	sHashes ( argv[1] , n , hash , preN , lastN , NULL , NULL , 0 ) ; 
	
	int i , j ; 
	int *idx = (int*) malloc( (N-k+1) * sizeof(int) ) ; 
	for( i = 0 ; i < N-k+1 ; i++ ) 
		idx[i] = n*i ; 
	
	sHashQuickSort ( hash , idx , n , N-k+1 ) ; 
	
	int N1 = strlen(argv[2]) ; 
	
	ull *hash2 = (ull*) malloc( n * (N1-k+1) * sizeof(ull) ) ; 
	int hit = sHashes ( argv[2] , n , hash2 , preN , lastN , hash , idx , N-k+1 ) ; 
	if( hit > 0 ) 
		printf( "a hit was found\n" ) ; 
	else 
		printf( "a hit was not found\n" ) ; 
	
	/* 
	for( i = 0 ; i < N-k+1 ; i++ ) 
	{
		for( j = 0 ; j < n ; j++ ) 
			printf( "%llu " , hash[ j + idx[i] ] ) ; 
		printf( "\n" ) ; 
	}
	*/ 
	
	free( hash ) ; 
	
	/*
	ull *hash = (ull*) malloc( n * sizeof(ull) ) ; 
	sHash( argv[1] , n , m , hash ) ; 
	printf( "%s\nHas hash\n" , argv[1] ) ; 
	int i ; 
	for( i = 0 ; i < n ; i++ ) 
		printf( "%llu" , hash[i] ) ; 
	printf( "\n" ) ; 
	*/
	
	return 0 ; 
}

