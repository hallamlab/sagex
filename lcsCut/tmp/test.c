
#include <stdio.h> 
#include <stdlib.h> 

#include "lcs.h"

int main() 
{
	char a1[] = "sdfkvbfgh" ; 
	char a2[] = "trjkvu1234567890bhkjfgh" ; 
	char a3[] = "hfdjaehsfsafdsag" ; 
	char b1[] = "qwgdsasfdn" ; 
	char b2[] = "xcsdfwe" ; 
	char b3[] = "cvnmfghjkwertyuio" ; 
	char b4[] = "ruiod1234567890fjkvbnm" ; 
	char b5[] = "wertysdfjkcvbnm" ; 
	char *a[3] ; 
	char *b[5] ; 
	a[0] = a1 ; 
	a[1] = a2 ; 
	a[2] = a3 ; 
	b[0] = b1 ; 
	b[1] = b2 ; 
	b[2] = b3 ; 
	b[3] = b4 ; 
	b[4] = b5 ; 
	
	int subset[5] ; 
	int subN ; 
	lcsSubset ( a , 3 , b , 5 , 4 , subset , &subN ) ; 
	
	printf( "%i\n" , subset[0] ) ; 
	
	return 0 ; 
}

