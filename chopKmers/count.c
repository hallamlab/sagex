
#include "count.h"

// Special thanks to Kishori M. Konwar for designing this function 

// calculates kmers in O(n) time 
// nucs : a string, should only contain the characters 'A', 'T', 'C', & 'G' 
// maxLen : if > 0 , only the first maxLen characters of nucs will be traversed, otherwise all.   
// out : a pre-allocated int array of 2^k entries
void kmerCounter ( char *nucs , int maxLen , int *out ) 
{
	int s ; 
    	for(int s=0; s < 256; s++)
		out[s]=0;
    	
    	unsigned int indx = 0, j =0; 
	
	int i ; 
        for( i = 0 ; *nucs != '\0' && ( maxLen < 1 || i < maxLen ) ; i++ )  
	{
       		switch(*nucs) 
		{
           		case 'A':
           		case 'a':
                      		indx = 4*indx + 0;  
                      		break;
           		case 'T':
           		case 't':
                      		indx = 4*indx + 1;  
                      		break;
           		case 'C':
           		case 'c':
                      		indx = 4*indx + 2;  
                      		break;
           		case 'G':
           		case 'g':
                      		indx = 4*indx + 3;  
                      		break;
			case 'N': 
			case 'n': 
				break ; // skip unknowns 
           		default:
                      	j=0;
			fprintf( stderr , "ERROR: kmerCounter: invalid character: %c, %s\n" , *nucs , nucs ) ; 
			return ; 
       		}   
       		j++;
       		nucs++;
       		indx = indx % 256;
       		if( j > 3)
			out[indx]++;
        } // end while
}

