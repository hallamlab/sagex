
#include "lcs.h"

// finds the length of the longest common substring 
// a : a character array of length an 
// b : a character array of length bn 
// tmp1 , tmp2 : workingspace of length max(an,bn) 
int dpLcs ( char *a , int an , char *b , int bn , int *tmp1 , int *tmp2 ) 
{
	int *prev = tmp1 ; 
	int *post = tmp2 ; 
	
	int max = 0 ; 
	int i , j ; 
	for( j = 0 ; j < bn ; j++ ) 
		prev[j] = 0 ; 
	for( i = 0 ; i < an ; i++ ) 
	{
		for( j = 0 ; j < bn ; j++ ) 
		{
			if( a[i] == b[j] ) 
			{
				if( i == 0 || j == 0 ) 
					post[j] = 1 ; 
				else 
					post[j] = prev[j-1] + 1 ; 
				
				if( post[j] > max ) 
					max = post[j] ; 
			}
			else 
				post[j] = 0 ;  
		} 
		tmp1 = post ; 
		post = prev ; 
		prev = tmp1 ; 
	} 
	return max ; 
}

