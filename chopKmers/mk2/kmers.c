
#include "count.h" 

struct kmerArg 
{
	
}; 

// POSIX thread wrapper 
void *pCount( void *arg ) 
{
	
}

void countKmers ( char **fasta , int n , int minLength , int chopSize , int overlap ) 
{
	// check assumptions 
	if( minLength < chopSize ) 
	{
		fprintf( stderr , "ERROR: minimum contig length must be at least as large as chop size\n" ) ; 
		return ; 
	}
	if( overlap >= chopSize  ) 
	{
		fprintf( stderr , "ERROR: overlap must be less than chopSize\n" ) ; 
		return ; 
	}
	
	// record contig lengths 
	int *lengths = (int*) malloc( n * sizeof(int) ) ; 
	int i ; 
	for( i = 0 ; i < n ; i++ ) 
		lengths[i] = strlen( fasta[i] ) ; 
	
	// extract subset of sufficiently long contigs 
	int *sub = (int*) malloc( n * sizeof(int) ) ; 
	int subN = -1 ; 
	subMinLen ( lengths , n , minLength , sub , &subN ) ; 
	
	// calculate total work to be done 
	int workN = -1 ; 
	totalWorkItems( sub , subN , lengths , chopSize , overlap , &workN ) ; 
	
	// order work 
	
	// record work items 
	
	// divide tasks 
	
	// populate arguments 
	
	// run threads 
	
	// join threads 
	
	// free memory 
	free( lengths ) ; 
	free( sub ) ; 
}

