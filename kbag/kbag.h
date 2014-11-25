
#include <limits.h> 
#include <math.h> 
#include <stdio.h> 
#include <stdlib.h> 

typedef unsigned long long ull ; 

// quick sorts a list of sHashes 
// list : storage space for the ulls 
// idx : storage place for the ull IDs 
// listN : number of sHashes in list 
// n : number of ull's per sHash 
void sHashQuickSort ( ull *list , int *idx , int n , int listN ) ; 

// get number of ull's required to store a key for a kmer of size k  
// k : kmer length 
// n : output , number of ull's required to store a key for a kmer of size k 
// preN : number of characters to represent in the first n-1 hash chunks 
// lastN : number of character to represent in the n-th hash chunk    
void sHashGetKeySize ( int k , int *n , int *preN , int *lastN ) ; 

// calculates the hashes for all substrings of length cut 
// a string of length N will have N - cut hashes 
// s : the string to be have substrings hashed 
// n : the number of ull's required to store a hash of size k  
// dict : if not null, sHashes will perform lookups.  
// sHashes will return 1 if a kmer from s is in dict , 0 otherwise. 
// idx : must not be null if dict is not null. Stores the sorted indices of dict 
// idxN : length of idx 
// sHashes returns -1 on errors  
int sHashes ( char *s , int n , ull *hash , int preN , int lastN , ull *dict , int *idx , int idxN ) ; 


