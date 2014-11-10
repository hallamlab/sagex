
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <pthread.h> 

// calculates kmers in O(n) time 
// nucs : a string, should only contain the characters 'A', 'T', 'C', & 'G' 
// maxLen : if > 0 , only the first maxLen characters of nucs will be traversed, otherwise all.   
// out : a pre-allocated int array of 2^k entries
void kmerCounter ( char *nucs , int maxLen , int *out ) ; 

// Subsets contigs by a minimum length parameter 
// legnths : n-length vector of contig lengths 
// cut : the minimum length for sequences 
// out : output, space for n integers to store the indecies of the subset 
// outN : output, number of indices in the subset 
void subMinLen ( int *lengths , int n , int cut , int *out , int *outN ) ; 

// calculates the total number of work items 
// sub : n-length array, a subset of contig indices representing sufficiently long contigs 
// lengths : lengths of all contigs 
// out : output, total work to be done 
void totalWorkItems( int *sub , int n , int *lengths , int chopSize , int overlap , int *out ) ; 

// records cumulative work per contig  
// use when not chopping contigs 
// sub : n-sized subset of indices of counted contigs 
// lengths : lengths of all contigs 
// chopSize : if < 1 , work will be the entire contig length, else it will be of chopSize per task 
// out : output, n-length list of cumulative work 
void getCumulativeWork ( int *sub , int n , int *lengths , int chopSize , int *out ) ; 

// records works items for threads 
// behaviour depends on chopSize 
// if workN == n , each contig in the subset (sub) gets a single work item 
// otherwise, each contig has at least one work item 
// sub : n-length list of indices of contig subset 
// lengths : lengths of all contigs 
// chopSize : must be greater than overlap 
// seq : output, workN-length array of task sequences 
// loc : output, workN-length array of task sequence locations  
void recordWorkItems( int *sub , int n , int workN , int *lengths , int chopSize , int overlap , int* seq , int* loc ) ; 

// for dividing tasks amongst threads 
// work : workN-length array of work accumulation 
// starts : output, threads-length array to store initial work items, inclusive 
// ends : output, threads-length array to store final work items, non-inclusive   
void divideTasks ( int *work , int workN , int *starts , int *ends , int threads ) ; 


