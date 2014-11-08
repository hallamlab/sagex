
#include "count.h"

// Subsets contigs by a minimum length parameter 
// legnths : n-length vector of contig lengths 
// cut : the minimum length for sequences 
// out : output, space for n integers to store the indecies of the subset 
// outN : output, number of indices in the subset 
void subMinLen ( int *lengths , int n , int cut , int *out , int *outN ) 
{
	int i ; 
	*outN = 0 ; 
	for( i = 0 ; i < n ; i++ ) 
	{
		if( lengths[i] >= cut ) 
		{
			out[ *outN ] = i ; 
			*outN += 1 ; 
		}
	}
}

