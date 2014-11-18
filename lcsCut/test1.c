
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
	
	int i, j ; 
	for( i = 0 ; i < 50 ; i++ ) 
	{
		a[i] = (char*) malloc( 100 * sizeof(char) ) ; 
		for( j = 0 ; j < 99 ; j++ ) 
			a[i][j] = 'a' ; 
		a[i][99] = '\0' ; 
	}
	for( i = 0 ; i < 200 ; i++ ) 
	{
		b[i] = (char*) malloc( 100 * sizeof(char) ) ; 
		for( j = 0 ; j < 99 ; j++ ) 
			b[i][j] = 'b' ; 
		b[i][99] = '\0' ; 
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
}

