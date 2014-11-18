
#include "lcs.h" 

struct lcsArg 
{
	char **ref ; 
	int refN ; 
	char **str ; 
	int start ; 
	int stop ; 
	int cut ; 
	int *subset ; 
	int *subN ; 
}; 

void *lcsFunc ( void *parg ) 
{
	struct lcsArg *arg = (struct lcsArg*) parg ; 
	lcsSubset ( arg->ref , arg->refN , &(arg->str[ arg->start ]) , arg->stop - arg->start , arg->cut , &( arg->subset[ arg->start ] ) , arg->subN ) ; 
	int i ; 
	for( i = 0 ; i < *(arg->subN) ; i++ ) 
		arg->subset[ arg->start + i ] += arg->start ; 
}

void lcsPosix ( char **ref , int refN , char **str , int strN , int cut , int* subset , int *subN , int threads ) 
{
	if( threads < 2 ) 
	{
		lcsSubset ( ref , refN , str , strN , cut , subset , subN ) ; 
		return ; 
	}
	
	int *work = (int*) malloc( strN * sizeof(int) ) ; 
	work[0] = strlen( str[0] ) ; 
	int i ; 
	for( i = 1 ; i < strN ; i++ ) 
		work[i] = work[i-1] + strlen( str[i] ) ; 
	int total = work[strN-1] ; 
	int div = total / threads ; 
	int prev = 0 ; 
	int *starts = (int*) malloc( threads * sizeof(int) ) ; 
	int *stops = (int*) malloc( threads * sizeof(int) ) ; 
	starts[0] = 0 ; 
	int thr = 0 ; 
	for( i = 0 ; i < strN ; i++ ) 
	{
		if( work[i] - prev > div ) 
		{
			stops[thr] = i+1 ; 
			prev = work[i] ; 
			thr++ ; 
			starts[thr] = i+1 ; 
		}
		if( thr + 1 == threads ) 
			break ; 
	}
	stops[thr] = strN ; 
	thr++ ; 
	threads = thr ; 
	
	pthread_t *pThreads = (pthread_t*) malloc( threads * sizeof(pthread_t) ) ; 
	struct lcsArg *parg = (struct lcsArg*) malloc( threads * sizeof(struct lcsArg) ) ; 
	int *N = (int*) malloc( threads * sizeof(int) ) ; 
	
	int pErr ; 
	for( i = 0 ; i < threads ; i++ ) 
	{
		parg[i].ref = ref ; 
		parg[i].refN = refN ; 
		parg[i].str = str ; 
		parg[i].start = starts[i] ; 
		parg[i].stop = stops[i] ; 
		parg[i].cut = cut ; 
		parg[i].subset = subset ; 
		parg[i].subN = &N[i] ; 
		
		pErr = pthread_create ( &pThreads[i] , NULL , lcsFunc , (void*) &parg[i] ) ; 
		if( pErr ) 
		{
			fprintf( stderr , "ERROR: longest common substring, posix\n" ) ; 
			return ; 
		}
	}
	
	void *status ; 
	for( i = 0 ; i < threads ; i++ ) 
	{
		pErr = pthread_join( pThreads[i] , &status ) ; 
		if( pErr ) 
		{
			fprintf( stderr , "ERROR: longest common substring, posix 2\n" ) ; 
			return ; 
		}
	}
	
	// reconstruct work 
	*subN = 0 ; 
	int j ; 
	int k = 0 ; 
	for( i = 0 ; i < threads ; i++ ) 
	{
		*subN += N[i] ; 
printf( "N[%i] = %i\n" , i , N[i] ) ; 
		for( j = 0 ; j < N[i] ; j++ ) 
		{
			subset[k] = subset[ starts[i] + j ] ; 
			k++ ; 
		}
	}
	
	free( work ) ; 
	free( starts ) ; 
	free( stops )  ; 
	free( pThreads ) ; 
	free( parg ) ; 
	free( N ) ; 
}











