
#include "count.h" 

struct kmerArg 
{
	int threadIdx ; 
	char **fasta ; 
	int chopSize ; 
	int *seq ; 
	int *loc ; 
	int *starts ; 
	int *ends ; 
	int *tmp ; // working space 
	int N ; 
	int *out ; // workN X 2^k column-major integer matrix 
}; 

// POSIX thread wrapper 
void *pCount( void *arg ) 
{
	struct kmerArg karg = (struct kmerArg) *arg ; 
	int i , j ; 
	for( i = karg.starts[ karg.threadIdx ] ; i < karg.ends[ karg.threadIdx ] ; i++ ) 
	{
		kmerCounter ( &(karg.fasta[ karg.seq[i] ][ karg.loc[i] ]) , karg.chopSize , karg.tmp ) ; 
		for( j = 0 ; j < 256 ; j++ ) 
			out[ i + N * j ] = karg.tmp[j] ; 
	}
}

void countKmers ( char **fasta , int n , int minLength , int chopSize , int overlap , int threads , int **out , int *outN ) 
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
	
	// allocate output space 
	*outN = workN ; 
	*out = (int*) malloc( workN * 256 * sizeof(int) ) ; 
	
	// cumulatively total work  
	int *work = (int) malloc( subN * sizeof(int) ) ; 
	getCumulativeWork ( sub , subN , lengths , chopSize , work ) ; 
	
	// record work items 
	int *seq = (int*) malloc( workN * sizeof(int) ) ; 
	int *loc = (int*) malloc( workN * sizeof(int) ) ; 
	recordWorkItems( sub , subN , workN , lengths , chopSize , overlap , seq , loc ) ; 
	
	// divide tasks 
	int *starts = (int*) malloc( threads * sizeof(int) ) ; 
	int *ends = (int*) malloc( threads * sizeof(int) ) ; 
	divideTasks ( work , workN , starts , ends , threads ) ; 
	
	// populate arguments 
	struct kmerArg *karg = (struct kmerArg*) malloc( threads * sizeof(struct kmerArg) ) ; 
	int **tmp = (int**) malloc( threads * sizeof(int*) ) ; 
	for( i = 0 ; i < threads ; i++ ) 
	{
		tmp[i] = (int*) malloc( 256 * sizeof(int) ) ; 
		karg[i].threadId = i ; 
		karg[i].fasta = fasta ; 
		karg[i].chopSize = chopSize ; 
		karg[i].seq = seq ; 
		karg[i].loc = loc ; 
		karg[i].starts = starts ; 
		karg[i].ends = ends ; 
		karg[i].tmp = tmp[i] ; 
		karg[i].N = workN ; 
		karg[i].out = out ; 
	}
	
	// run threads 
	pthread_t *pThreads = (pthread_t*) malloc( threads * sizeof(pthread_t) ) ; 
	int pErr ; 
	for( i = 0 ; i < threads ; i++ ) 
	{
		pErr = pthread_create ( &pThreads[i] , NULL , pCount , (void*) &karg[i] ) ; 
		if( pErr )
		{
			fprintf( stderr , "ERROR, POSIX: %i\n" , pErr ) ; 
			return ;
		}
	}
	
	// join threads 
	void *status ; 
	for( i = 0 ; i < threads ; i++ ) 
	{
		pErr = pthread_join ( pThreads[i] , &status ) ; 
		if( pErr ) 
		{
			fprintf( stderr , "ERROR: POSIX: %i\n" , pErr ) ; 
			return ; 
		}
		free( tmp[i] ) ;  
	}
	
	// free memory 
	free( lengths ) ; 
	free( sub ) ; 
	free( work ) ; 
	free( seq ) ; 
	free( loc ) ; 
	free( kmerArg ) ; 
	free( starts ) ; 
	free( ends ) ; 
	free( tmp ) ; 
}




