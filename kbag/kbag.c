
#include <stdio.h> 
#include <stdlib.h> 
#include <limits.h> 
#include <string.h> 

typedef unsigned long long ull ; 

struct kbag 
{
	ull n ; 
	ull logn ; // for O( log log n ) lookups 
	ull *bag ; 
}; 

void counter ( char *nucs , int *out , int k )  
{
        ull indx = 0, j = 0 ; 
    
        int i ; 
        for( i = 0 ; *nucs != '\0' ; i++ )   
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

void makeBag ( char **fasta , int n , kbag **bag , int k ) 
{
	int *lengths = (int*) malloc( n * sizeof(int) ) ; 
	int i ; 
	ull total = 0 ; 
	for( i = 0 ; i < n ; i++ ) 
	{ 
		lengths[i] = strlen( fasta[i] ) ; 
		if( lengths[i] >= k ) 
			total += (ull) (lengths[i] - k-1) ; 
	} 
	if( total == 0 ) 
	{
		fprintf( stderr , "WARNING: makeBag: empty bag!\n" ) ; 
		return ; 
	}
	
	
	
	free( lengths ) ; 
}

