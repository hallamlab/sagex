
#include "chopCounter.h"

// counts occurrence of a substring in a super string 
int count ( char *sub , char *sup , int n )
{
        // printf( "DEBUG: comparing %s and %s\n" , sub , sup ) ; 
        if( strlen( sub ) == 0 )
                return -1 ;
        if( strlen( sup ) == 0 )
                return 0 ;
        if( strlen( sub ) > strlen( sup ) )
                return 0 ;
	
        int out = 0 ;
        int i , j , flag ;
        for( i = 0 ;  i < n ; i++ )
        {
                flag = 1 ;
                for( j = 0 ; flag > 0 && sub[j] != '\0' && sup[i+j] != '\0' ; j++ )
                {
                        if( sup[i+j] != sub[j] )
                                flag = 0 ;
                }
                if( sup[i+j] == '\0' && sub[j] != '\0' )
                        flag = 0 ;
                if( flag > 0 )
                        out = out + 1 ;
        }
        return out ;
}

