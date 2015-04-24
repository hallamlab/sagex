
#include "count.h" 

// assumes length >= chopSize 
int countChops ( int length , int chopSize , int overlap ) 
{
	if( chopSize > length )
		return 0 ; 
        int m = chopSize - overlap ; 
	int out = floor( ((double) (length - chopSize) ) / ((double) m) ) + 1 ; 
	return out ; 
}

// calculates the total number of work items 
// sub : n-length array, a subset of contig indices representing sufficiently long contigs 
// lengths : lengths of all contigs 
// out : output, total work to be done 
void totalWorkItems( int *sub , int n , int *lengths , int chopSize , int overlap , int *out ) 
{
	if( chopSize < 1 ) 
	{
		*out = n ; 
		return ; 
	}
	*out = 0 ; 
	int i ; 
	for( i = 0 ; i < n ; i++ ) 
		*out += countChops( lengths[ sub[i] ] , chopSize , overlap ) ; 
}

// records cumulative work per contig  
// use when not chopping contigs 
// sub : n-sized subset of indices of counted contigs 
// lengths : lengths of all contigs 
// out : output, n-length list of cumulative work 
void getCumulativeWork ( int *sub , int n , int *lengths , int *out ) 
{
	int i ; 
	out[0] = lengths[ sub[0] ] ; 
	for( i = 1 ; i < n ; i++ ) 
		out[i] = out[i-1] + lengths[ sub[i] ] ; 
	return ; 
}

// records works items for threads 
// behaviour depends on chopSize 
// if workN == n , each contig in the subset (sub) gets a single work item 
// otherwise, each contig has at least one work item 
// sub : n-length list of indices of contig subset 
// lengths : lengths of all contigs 
// chopSize : must be greater than overlap 
// seq : output, workN-length array of task sequences 
// loc : output, workN-length array of task sequence locations  
void recordWorkItems( int *sub , int n , int workN , int *lengths , int chopSize , int overlap , int* seq , int* loc ) 
{
	int i , j , k ; 
	if( workN == n ) // not chopping contigs  
	{
		for( i = 0 ; i < n ; i++ ) 
		{
			seq[i] = sub[i] ; 
			loc[i] = 0 ; // sub[i] ; 
		}
		return ; 
	}
	j = 0 ; // work item number 
	k = 0 ; // location on contig 
	for( i = 0 ; i < n ; i++ ) // contig number 
	{
		// record items for sequence i 
		while( k + chopSize <= lengths[ sub[i] ] ) 
		{
			if( j >= workN ) 
				fprintf( stderr , "ERROR: chopKmers: recordWorkItems: Allocation overflow!\n" ) ; 
			
			seq[j] = sub[i] ; // NOTE THAT SUB IS NO LONGER NEEDED AFTER THIS POINT 
			loc[j] = k ; 
			j++ ; 
			
			k +=  chopSize - overlap ; 
		}
		
		// increment sequence 
		k = 0 ; 
	}
	if( j != workN ) 
		fprintf( stderr , "WARNING: chopKmers: recordWorkItems: Allocation incomplete!\n" ) ; 
}

// for dividing tasks amongst threads 
// work : workN-length array of work accumulation 
// chopSize 
// starts : output, threads-length array to store initial work items, inclusive 
// ends : output, threads-length array to store final work items, non-inclusive 
// returns actual number of threads needed  
int divideTasks ( int *work , int workN , int chopSize , int *starts , int *ends , int threads ) 
{
	int i ; 
	if( workN == threads ) // easy case  
	{
		starts[0] = 0 ; 
		for( i = 0 ; i < workN-1 ; i++ ) 
		{
			starts[i] = i ; 
			ends[i] = i + 1 ; 
		}
		ends[workN-1] = workN ; 
		return threads ; 
	}
	int div ; 
	if( chopSize > 0 ) // if chopping contigs 
	{
		if( workN % threads == 0 ) 
			div = workN / threads ; 
		else 
			div = (workN + 1) / threads ; 
		for( i = 0 ; i < threads ; i++ ) 
		{
			starts[i] = i*div ; 
			ends[i] = (i+1)*div ; 
		}
		ends[ threads-1 ] = workN ;  
		return threads ; 
	}
	// if not chopping contigs 
	if( work[ workN-1 ] % threads == 0 ) 
		div = work[ workN-1 ] / threads ; 
	else 
		div = (work[ workN-1 ] + 1) / threads ; 
	starts[0] = 0 ; 
	int prev = 0 ; 
	int thr = 0 ; 
	for( i = 0 ; i < workN ; i++ ) 
	{
		if( work[i] - prev >= div ) 
		{
			prev = work[i] ; 
			if( thr + 1 == threads ) 
			{
				ends[ thr ] = workN ; 
				return thr+1 ; 
			}
			else 
			{
				ends[ thr ] = i ; 
				thr += 1 ; 
				starts[ thr ] = i ; 
			}
		}
	}
	ends[thr] = workN ; 
	return thr+1 ; 
}











