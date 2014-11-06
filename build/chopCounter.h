
#include <pthread.h>

#include <string.h> 
#include <stdio.h> 
#include <stdlib.h>

int count ( char *sub , char *sup , int n ) ; 

void tetraCounterChop ( char *nucs , int *out , int chopSize ) ; 

// fasta : an array of strings, only ATCG characters allowed  
// n : the number of strings 
// minLen : minimum contig length  
// out : a non-allocated pointer to an int pointer to store a 256 X ? matrix  
// outN : the number of columns of out 
// names : an output for enumerating the names of sequences, if NULL nothing will be output 
void posixCounter( char **fasta , int n , int minLen , int threads , int chopSize , int overlap , int **out , int *outN , int **names ) ; 

