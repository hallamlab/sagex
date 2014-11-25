
#include "lcs.h" 

int main( int argc , char **argv )
{
	if( argc < 2 ) 
	{
		printf( "provide number of threads\n" ) ; 
		return -1 ; 
	}
	int thr = atoi( argv[1] ) ; 
	
	char **a = (char**) malloc( 50 * sizeof(char*) ) ; 
	char **b = (char**) malloc( 200 * sizeof(char*) ) ; 
	
	char numstr[10] ; 
	int i, j ; 
	for( i = 0 ; i < 50 ; i++ ) 
	{
		sprintf( numstr , "%d" , i ) ; 
		a[i] = (char*) malloc( 110 * sizeof(char) ) ; 
		for( j = 0 ; j < 10 ; j++ ) 
			a[i][j] = 'a' ; 
		a[i][10] = '\0' ; 
		strcat( a[i] , numstr ) ;  
	}
	for( i = 0 ; i < 200 ; i++ ) 
	{
		sprintf( numstr , "%d" , i ) ; 
		b[i] = (char*) malloc( 110 * sizeof(char) ) ; 
		for( j = 0 ; j < 10 ; j++ ) 
			b[i][j] = 'b' ; 
		b[i][10] = '\0' ; 
		strcat( b[i] , numstr ) ; 
	}
	
	strcpy( a[5] , "xcdfhjqwuioqwertyuiopsdfghjvbn" ) ; 
	strcpy( b[143] , "xcvbnmdfghjertuizsdfgqwertyuiopsdxcvvb" ) ; 
	strcpy( b[122] , "cvnmfghjkqwertyuioghjvbndfgqwertyuio" ) ; 
	
	int subN ; 
	int *subset = (int*) malloc( 200 * sizeof(int) ) ; 
	
	lcsPosix ( a , 50 , b , 200 , 4 , subset , &subN , thr ) ; 
	
	for( i = 0 ; i < subN ; i++ )
	{
		printf( "%i\n" , subset[i] ) ; 
	}
	
	for( i = 0 ; i < 50 ; i++ ) 
		printf( "%i: %s\n" , i , a[i] ) ; 
	for( i = 0 ; i < 200 ; i++ ) 
		printf( "%i: %s\n", i , b[i] ) ; 
	
}



