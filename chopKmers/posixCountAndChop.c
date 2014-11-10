#include <math.h>
#include "chopCounter.h"

struct parg // posix arg 
{
	char **fasta ; 
	int start ; // index of first working object 
	int end ; // index immediately following last working object 
	int **out ; // the output matrix, to be allocated by the thread 
	int *cols ; // the thread will calculate the total number of chops 
	int chopSize ; 
	int overlap ;  
	int **names ; 
	int threadNum ; 
	int minLen ;  
}; 

int countChops ( char *str , int chopSize , int overlap ) 
{
	int m = chopSize - overlap ; 
	int out = floor( ((double) (strlen(str) - chopSize) ) / ((double) m) ) ; 
	if( out < 0 ) 
		return 0 ; 
	else
		return out ; 
}

void *pCount ( void *targ ) // temp arg  
{	
	struct parg *arg = (struct parg*) targ ; 
	int i , tmp ; 
	
	// count the total number of chops and therefore the number of columns to allocate 
	int n = 0 ; 
	if( arg->chopSize > 0 ) // number of kmer vectors will be greater than the number of contigs  
	{
		for( i = arg->start ; i < arg->end ; i++ ) 
		{
			if( strlen( arg->fasta[i] ) > arg->minLen ) 
			{
				tmp = countChops( arg->fasta[i] , arg->chopSize , arg->overlap ) ; 
				n += tmp ; 
			}
		}
	}
	else // at most one kmer per contig 
	{
		for( i = arg->start ; i < arg->end ; i++ ) 
		{
			if( strlen( arg->fasta[i] ) > arg->minLen ) 
				n++ ; 
		}
	}
	
	*(arg->cols) = n ; 
	
	// allocate the space 
	*(arg->out) = (int*) malloc( 256 * n * sizeof(int) ) ; 
	
	int m = arg->chopSize - arg->overlap ; // adjustment size 
	
	// allocate space for names, if requested 
	if( arg->names != NULL ) 
		*(arg->names) = (int*) malloc( n * sizeof(int) ) ; 
	
	i = 0 ; // chop number overall 
	int j = arg->start ; // sequence number 
	int k = 0 ; // chop number in current sequence 
	tmp = strlen( arg->fasta[ j ] ) ; 
	if( arg->chopSize > 0 ) // chop up contigs  
	{ 
		// skip all short contigs  
		while( strlen( arg->fasta[i] ) > arg->minLen && i < n ) 
			i++ ; 
		while( i < n )  
		{
			// skip all short contigs 
			while( strlen( arg->fasta[i] ) > arg->minLen && i < n ) 
				i++ ; 
			
			if( k * m + arg->chopSize > tmp ) // switch to next sequence 
			{
				j = j + 1 ; 
				k = 0 ; 
				tmp = strlen( arg->fasta[ j ] ) ; 
			}
			
			tetraCounterChop( &(arg->fasta[j][k*m]) , &( (*(arg->out))[ 256 * i ] ) , arg->chopSize ) ; 
// int l ; int flag = 0 ; for( l = 0 ; l < 256 && flag == 0 ; l++ ){ if( (*(arg->out))[ l + 256 * i ] > 0 ){ flag = 1 ; } } ; if( flag == 0 ) fprintf( stderr , "DEBUG!\n" ) ; 
			
			// name the chop, if requested 
			if( arg->names != NULL ) 
				(*arg->names)[i] = j ; 
			
			// iterate chop number within current sequence 
			k++ ; 
			
			// iterate overal chop number 
			i++ ; 
		}
	}
	else // do not chop up contigs 
	{ 
		j = 0 ; 
		for( i = arg->start ; i < arg->end ; i++ ) // one kmer vector per contig  
		{
			if( strlen( arg->fasta[i] ) > arg->minLen ) 
			{
				tetraCounterChop( arg->fasta[i] , &( (*(arg->out))[ 256 * j ] ) , strlen( arg->fasta[i] ) ) ; 
				
				// name the contig, if requested 
				if( arg->names != NULL ) 
					(*arg->names)[j] = i ; 
				j++ ; 
			} 
		}
	}
	pthread_exit(NULL) ; 
}

// fasta : an array of strings, only ATCG characters allowed  
// n : the number of strings 
// minLen : minimum contig length  
// out : a non-allocated pointer to an int pointer to store a 256 X ? matrix  
// outN : the number of columns of out 
// names : an output for enumerating the names of sequences, if NULL nothing will be output 
void posixCounter( char **fasta , int n , int minLen , int threads , int chopSize , int overlap , int **out , int *outN , int **names ) 
{
	if( threads > n )
		threads = n ; 
	pthread_t *thr = (pthread_t*) malloc( threads * sizeof(pthread_t) ) ; 
	struct parg *pargs = (struct parg*) malloc( threads * sizeof(struct parg) ) ; 
	
	int **mats = (int**) malloc( threads * sizeof(int*) ) ; 
	int *cols = (int*) malloc( threads * sizeof(int) ) ; 
	int **nameTemp = NULL ; 
	if( names != NULL ) 
		nameTemp = (int**) malloc( threads * sizeof(int*) ) ; 
	
	// design work 
	int i ; 
	for( i = 0 ; i < threads ; i++ ) 
	{
		pargs[i].fasta = fasta ; 
		pargs[i].start = i * n / threads ; // TODO this is a terrible way to distribute work   
		pargs[i].end = (i+1) * n / threads ; 
		if( i + 1 == threads ) 
			pargs[i].end = n ; 
		pargs[i].chopSize = chopSize ;  
		pargs[i].overlap = overlap ; 
		pargs[i].out = &mats[i] ; 
		pargs[i].cols = &cols[i] ; 
		pargs[i].threadNum = i ; 
		pargs[i].minLen = minLen ; 
		
		if( names != NULL ) 
			pargs[i].names = &nameTemp[i] ; 
		else
			pargs[i].names = NULL ; 
	}
	
	// start work 
	int err ; 
	for( i = 0 ; i < threads ; i++ ) 
	{
		// printf( "DEBUG: Starting thread %i, with start: %i, & end: %i\n" , i , pargs[i].start , pargs[i].end ) ; 
		
		err = pthread_create( &thr[i] , NULL , pCount , (void*) &pargs[i] ) ; 
		if( err ) 
		{
			fprintf( stderr , "ERROR, POSIX: %i\n" , err ) ; 
			return ; 
		}
	}
	
	// join threads 
	void *status ; 
	for( i = 0 ; i < threads ; i++ ) 
	{
		err = pthread_join( thr[i] , &status ) ; 
		if( err ) 
		{
			fprintf( stderr , "ERRPR, POSIX: %i\n" , err ) ; 
			return ; 
		}
	}
	
	// get the total number of columns 
	*outN = 0 ; 
	for( i = 0 ; i < threads ; i++ ) 
	{
		*outN += cols[i] ; 
	}
	// allocate the out-going matrix 
	*out = (int*) malloc( *outN * 256 * sizeof(int) ) ; 
fprintf( stderr , "%i %i %i %i\n" , (*out)[0], (*out)[1], (*out)[2], (*out)[3] ) ; // prints zeros !!!  
	
	// allocate name space, if requested 
	if( names != NULL ) 
		*names = (int*) malloc( *outN * sizeof(int) ) ; 
	
	// copy each submatrix into the out-going matrix 
	int j , k ; 
	int m = 0 ; 
	int adj = 0 ; 
	for( i = 0 ; i < threads ; i++ ) 
	{
		for( j = 0 ; j < cols[i] ; j++ ) // cycle thru the columns of the i-th kmer matrix made by the i-th thread 
		{
			for( k = 0 ; k < 256 ; k++ ) 
				(*out)[ m + *outN * k ] = mats[i][ k + 256 * j ] ; // notice transpose  
//for( k = 0 ; k < 256 ; k++ ) fprintf( stderr , "%i " , mats[i][ k + 256 * j ] ) ; fprintf( stderr , "\n" ) ; // zeros in mats[i][*] 
			
			// record the name, if requested 
			if( names != NULL ) 
				(*names)[ m ] = nameTemp[i][j] ; 
			m++ ; 
		}
		
		free( mats[i] ) ; 
	}
// for( i = 0 ; i < *outN ; i++ ){int flag = 0; for( j=0 ;j < 256 && flag == 0 ; j++ ) { if( (*out)[ j + 256 * i ] > 0 ) flag = 1 ; } ; if( flag == 0 ) fprintf( stderr , "D!\n" ) ; } 
	
	free( thr ) ; 
	free( pargs ) ; 
	free( mats ) ; 
	free( cols ) ; 
	if( names != NULL ) 
		free( nameTemp ) ; 
	
	return ; 
}















