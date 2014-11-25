
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <pthread.h> 

// finds the length of the longest common substring 
// a : a character array of length an 
// b : a character array of length bn 
// tmp1 , tmp2 : workingspace of length max(an,bn) 
// cut : if > 0 , process will return when a substring of size cut is found 
int dpLcs ( char *a , int an , char *b , int bn , int *tmp1 , int *tmp2 , int cut ) ; 

// subsets a set of strings (str) by minimum longest common substrings 
// ref : the reference genome of refN-strings. 
// str : the subject for subsetting 
// cut : the minimum length for a longest common substring 
// subset : output space, length strN, contains indices of the output subset of strings 
// subN : output, a single int, length of the subset 
void lcsSubset ( char **ref , int refN , char **str , int strN , int cut , int* subset , int *subN ) ; 

void lcsPosix ( char **ref , int refN , char **str , int strN , int cut , int* subset , int *subN , int threads ) ; 

