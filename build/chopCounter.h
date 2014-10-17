
#include <pthread.h>

#include <string.h> 
#include <stdio.h> 
#include <stdlib.h>

int count ( char *sub , char *sup , int n ) ; 

void tetraCounterChop ( char *nucs , int *out , int chopSize ) ; 

void posixCounter( char **fasta , int n , int threads , int chopSize , int overlap , int **out , int *outN , int **names ) ; 

