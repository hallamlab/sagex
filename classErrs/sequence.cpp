
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
using namespace std ; 

int main ( int argc , char** argv )
{
	if( argc < 4 )
	{
	    cout << "Simulates sequencing from a single chromosome FASTA. Outputs to stdout" << endl ; 
		cout << "Randomly assigns material to (1) stdout and (2) stderr" << endl ; 
	    cout <<  "Please provide (1) a fasta file , (2) the sequence size, and (3) the number of samples" << endl ; 
	    return -1 ; 
	}
	
	///////////////// LOAD IN THE SEQUENCE 
	
	int seqN = atoi( argv[2] ) ; 
	int sampN = atoi( argv[3] ) ; 
	
	char *seq = NULL ; 
	string s ; 
	ifstream fasta( argv[1] ) ; 
	
	getline( fasta , s ) ; 
	while( s[0] != '>' && fasta.good() )
		getline( fasta , s ) ; 
	getline( fasta , s ) ; 
	seq = (char*) malloc( ( s.length() + 1 ) * sizeof(char) ) ; 
	strcpy( seq , s.c_str() ) ; 
	
	while( fasta.good() ) 
	{
		int len = strlen( seq ) ; 
		getline( fasta , s ) ; 
		seq = (char*) realloc( seq , (len + s.length() + 1 ) * sizeof(char) ) ; 
		strcat( seq , s.c_str() ) ; 
	}
	fasta.close() ; 
	
	/////////////////// SAMPLE FROM THE SEQUENCE 
	
	int i , j ; 
	int N = strlen(seq) ; 
	char *buff = (char*) malloc( (seqN + 1) * sizeof(char) ) ; 
	for( i = 0 ; i < sampN ; i++ ) 
	{
		// stdout
		cout << ">Samp_" << i << endl ; 
		j = rand() % ( N - seqN ) ; 
		memcpy( buff , &seq[j] , seqN ) ; 
		buff[seqN] = '\0' ; 
		cout << buff << endl ; 
		
		// stderr 
		cerr << ">Samp_" << i << endl ; 
		j = rand() % ( N -  seqN ) ; 
		memcpy( buff , &seq[j] , seqN ) ; 
		buff[seqN] = '\0' ; 
		cerr << buff << endl ; 
	}
	
	return 1 ; 
}







