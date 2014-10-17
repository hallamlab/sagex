
#include "chopCounter.h"

// vec : a string, should only contain the characters 'A', 'T', 'C', & 'G' 
// out : a pre-allocated int array of 256 entries 
// chopSize : length of substring of nucs used to count kemrs 
void tetraCounterChop ( char *nucs , int *out , int chopSize ) 
{
        // printf( "DEBUG: counting for %s\n" , nucs ) ; 

        char gen[] = {'A' , 'T' , 'C' , 'G' } ;
        char tetra[]= { 'A' , 'A' , 'A' , 'A' , '\0' } ;
        int idx[] = { 0 , 0 , 0 , 0 } ;

        int i = 0 ;
        int j ;
        while( i < 256 )
        {
                // count tetra occurrences in string i 
                out[i] = count( &tetra[0] , nucs , chopSize ) ;

                // printf( "%s: %i\n" , tetra , out[i] ) ; // DEBUG 

                // iterate to next tetra 
                i = i + 1 ;
                idx[0] = idx[0] + 1 ;
                for( j = 0 ; j < 3 ; j++ )
                {
                        if( idx[j] == 4 )
                        {
                                idx[j] = 0 ;
                                idx[j+1] = idx[j+1] + 1 ;
                        }
                }
                for( j = 0 ; j < 4 ; j++ )
                        tetra[j] = gen[ idx[j] ] ;
        }
}




