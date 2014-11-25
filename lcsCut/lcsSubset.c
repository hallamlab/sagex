
#include "lcs.h" 

// subsets a set of strings (str) by minimum longest common substrings 
// ref : the reference genome of refN-strings. 
// str : the subject for subsetting 
// cut : the minimum length for a longest common substring 
// subset : output space, length strN, contains indices of the output subset of strings 
// subN : output, a single int, length of the subset 
void lcsSubset ( char **ref , int refN , char **str , int strN , int cut , int* subset , int *subN ) 
{
	int *refLengths = (int*) malloc( refN * sizeof(int) ) ; 
	int *strLengths = (int*) malloc( strN * sizeof(int) ) ; 
	
	int maxL = 0 ; 
	int i , j ; 
	for( i = 0 ; i < refN ; i++ ) 
	{
		refLengths[i] = strlen( ref[i] ) ; 
		if( refLengths[i] > maxL ) 
			maxL = refLengths[i] ; 
	}
	for( i = 0 ; i < strN ; i++ ) 
	{
		strLengths[i] = strlen( str[i] ) ; 
		if( strLengths[i] > maxL ) 
			maxL = strLengths[i] ; 
	}
	
	int *tmp1 = (int*) malloc( maxL * sizeof(int) ) ; 
	int *tmp2 = (int*) malloc( maxL * sizeof(int) ) ; 
	
	int flag ; 
	int len ; 
	*subN = 0 ; 
	for( i = 0 ; i < strN ; i++ ) 
	{
		flag = 0 ; 
		for( j = 0 ; j < refN && flag == 0 ; j++ ) 
		{
			len = dpLcs ( str[i] , strLengths[i] , ref[j] , refLengths[j] , tmp1 , tmp2 , cut ) ; 
			if( len >= cut ) 
			{
				flag = 1 ; 
				subset[ *subN ] = i ; 
				*subN += 1 ; 
			}
		}
	}
	
	free( refLengths ) ; 
	free( strLengths ) ; 
	free( tmp1 ) ; 
	free( tmp2 ) ; 
}

