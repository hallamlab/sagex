
#include "lcs.h"

// finds the length of the longest common substring 
// a : a character array of length an 
// b : a character array of length bn 
// tmp1 , tmp2 : workingspace of length max(an,bn) 
// if cut > 0 , process will cut out when a longest substring of size cut is found 
int dpLcs ( char *a , int an , char *b , int bn , int *tmp1 , int *tmp2 , int cut ) 
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
// fprintf( stderr , "i: %i, j: %i, an: %i, bn: %i ... " , i , j , an , bn ) ; 
// fprintf( stderr , "a[%i]: %i ," , i , a[i] ) ; 
// fprintf( stderr , " b[%i]: %i\n" , j , b[j] ) ; 
			if( a[i] == b[j] ) 
			{
				if( i == 0 || j == 0 ) 
					post[j] = 1 ; 
				else 
					post[j] = prev[j-1] + 1 ; 
				
				if( post[j] > max ) 
				{
					max = post[j] ; 
					if( cut > 0 && max >= cut ) 
						return max ; 
				}
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

