
#include <pthread.h>
#include <string.h> 
#include <stdio.h> 
#include <stdlib.h>
#include <string>


int count ( char *sub , char *sup ) ; 

int count ( char *sub , char *sup , int n ) ; 

bool kmerCounter ( char *nucs , unsigned int *out,  unsigned int k=4 ) ; 


// fasta : an array of strings, only ATCG characters allowed  
// n : the number of strings 
// out : a pre-allocated 256 X n column-indexed int matrix in array format 
void posixCounter ( char **fasta , int n , unsigned int *out , int threads ) ; 

void computeKmerVectors(std::string &content, unsigned int size, std::ostream *output);
void printVector(unsigned int *freq);
