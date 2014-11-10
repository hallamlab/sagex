
#include "utilities.h"
#include <string>
#include <iostream>
#include <exception>

#include <vector>
#define VSIZE 256
using namespace std;

// Generate tetra strings, search for them in a given block of strings 

// vec : a string, should only contain the characters 'A', 'T', 'C', & 'G' 
// out : a pre-allocated int array of 256 entries 
bool kmerCounter ( char *nucs , unsigned int *out, unsigned int k )  {

   // printf("<<<<\n");
    for(int s=0; s < 256; s++) out[s]=0;
    
    unsigned int indx = 0, j =0;
    bool amb = false;

	while( *nucs != '\0') {
       switch(*nucs) {
           case 'A':
           case 'a':
                      indx = 4*indx + 0; 
                      break;
           case 'T':
           case 't':
                      indx = 4*indx + 1; 
                      break;
           case 'C':
           case 'c':
                      indx = 4*indx + 2; 
                      break;
           case 'G':
           case 'g':
                      indx = 4*indx + 3; 
                      break;
           default:
                      j=0;
    //                  printf("Error  <%c>\n", *nucs);
                      amb = true;
                   //  printf("%d %d\n",j, indx);
       }  
       j++;
       nucs++;
       indx = indx % 256;
       if( j > 3) out[indx]++;
	} // end while
  //  printf(">>>>>>\n");
    return true;

}


void printVector(std::vector< std::vector<unsigned int> > vectors, std::ostream *output) {
     unsigned int i =0;
     for(std::vector< std::vector<unsigned int> >::iterator it= vectors.begin(); it != vectors.end(); it++) {
        i = 0;
        for(std::vector<unsigned int>::iterator it1= it->begin() ; it1 != it->end(); it1++, i++) {
           if(i!=0) *output << "\t";
           *output << *it1;
        }
        *output << "\n";
     }
}

void computeKmerVectors(std::string &content, unsigned int size, std::ostream *output) {
     static char *cstring = 0;
     unsigned int freq[VSIZE];
     static int i =0;
     std::vector< std::vector< unsigned int > > vectors;
     std::vector< unsigned int > simplevector;

     if(cstring==0) {
        cstring = (char *)malloc(1000000000*sizeof(char));
      //  printf("Allocating memory\n");

     }
    // printf("L=%d\n",strlen(content.c_str()));
     strcpy(cstring, content.c_str());
     char *start = cstring  + strlen(cstring) - size;

     while(cstring < start) { 
 
       //  printf("<%d\n",i);
         try {
            kmerCounter(start, freq);
         }
         catch(exception &e) {
            cout << "Error in kmer Counter " << e.what()  << std::endl;
         }

         simplevector.clear();
         for(int j=0; j < VSIZE; j++) simplevector.push_back(freq[j]);
         vectors.push_back(simplevector);

         if( vectors.size() > 1000 ) {
            printVector(vectors, output);
            vectors.clear();
         }

         if( ++i%1000 ==0) printf("%d\n",i);
         *start = '\0';
         start = start -  size;
     }

     if( vectors.size() > 0 ) {
        printVector(vectors, output);
        vectors.clear();
     }

}
// counts occurrence of a substring in a super string 
int count ( char *sub , char *sup ) 
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
	for( i = 0 ; sup[i+1] != '\0' ; i++ ) 
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



struct parg // posix arg 
{
	char **fasta ; 
	int start ; // index of first working object 
	int end ; // index immediately following last working object 
	unsigned int *out ; // the output matrix 
}; 

void *pCount ( void *targ ) // temp arg 
{
	struct parg *arg = (struct parg*) targ ; 
	int i ;
	// printf( "DEBUG: start %i, end: %i\n" , arg->start , arg->end ) ;  
	for( i = arg->start ; i < arg->end ; i++ ) 
	{
		kmerCounter( arg->fasta[i] , &( arg->out[ 256 * i ] ) ) ; 
	}
	pthread_exit(NULL) ; 
}

// fasta : an array of strings, only ATCG characters allowed  
// n : the number of strings 
// out : a pre-allocated 256 X n column-indexed int matrix in array format 
void posixCounter( char **fasta , int n , unsigned int *out , int threads ) 
{
	if( threads > n )
		threads = n ; 
	pthread_t *thr = (pthread_t*) malloc( threads * sizeof(pthread_t) ) ; 
	struct parg *pargs = (struct parg*) malloc( threads * sizeof(struct parg) ) ; 
	
	// design work 
	int i ; 
	for( i = 0 ; i < threads ; i++ ) 
	{
		pargs[i].fasta = fasta ; 
		pargs[i].start = i * n / threads ; 
		pargs[i].end = (i+1) * n / threads ; 
		if( i + 1 == threads ) 
			pargs[i].end = n ; 
		pargs[i].out = out ; 
		// printf( "DEBUG: Designing work for thread %i, start = %i, end = %i\n" , i , pargs[i].start , pargs[i].end ) ; 
	}
	
	// start work 
	int err ; 
	for( i = 0 ; i < threads ; i++ ) 
	{
		// printf( "DEBUG: Starting thread %i, with start: %i, & end: %i\n" , i , pargs[i].start , pargs[i].end ) ; 
		
		err = pthread_create( &thr[i] , NULL , pCount , (void*) &pargs[i] ) ; 
		if( err ) 
		{
			printf( "ERROR, POSIX: %i\n" , err ) ; 
			return ; 
		}
	}
	
	// join threads 
	void *status ; 
	for( i = 0 ; i < threads ; i++ ) 
	{
		err = pthread_join( thr[i] , &status ) ; 
		if( err ) 
		{
			printf( "ERRPR, POSIX: %i\n" , err ) ; 
			return ; 
		}
	}
	
	free( thr ) ; 
	free( pargs ) ; 
	
	return ; 
}

