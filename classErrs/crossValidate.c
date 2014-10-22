
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void crossValidate ( int *kk , int *pp , int *aa ) 
{
	int k = *kk ; 
	int p = *pp ; 
	int a = *aa ; 
	printf( "Simulating sequencing...\n" ) ; 
	system( "./sequence ./../../../mg1655eColi.fasta 4000 1000 1> tmp1 2> tmp2" ) ; // tmp1 : test , tmp2 : train 
	printf( "catting... \n" ) ; 
	system( "cat ./../../../metaG.fasta tmp1 > tmp3" ) ; // tmp3 : final test data set 
	printf( "making blast db...\n" ) ; 
	system( "makeblastdb -in tmp3 -dbtype nucl" ) ; 
	printf( "running blastn...\n" ) ; 
	system( "blastn -query ./tmp2 -db tmp3 -outfmt 6 > tmp4" ) ; // tmp4 : database 
	printf( "runnning sagex...\n" ) ; 
	char cmd[1000] ; 
	// cmd[0] = '\n' ; 
	strcpy( cmd , "../sagex -i tmp2 -G tmp3 -b tmp4 -t 50" ) ; 
	strcat( cmd , " -k " ) ; 
	char ks[100] ; 
	char ps[100] ; 
	char as[100] ; 
	sprintf( ks , "%i" , k ) ; 
	sprintf( ps , "%i" , p ) ; 
	sprintf( as , "%i" , a ) ; 
	strcat( cmd , ks ) ; 
	strcat( cmd , " -p " ) ; 
	strcat( cmd , ps ) ; 
	strcat( cmd , " -a " ) ; 
	strcat( cmd , as ) ; 
	strcat( cmd , " > tmp5" ) ; 
printf( "DEBUG, k: %i, p: %i, a: %i\n" , k , p , a ) ; 
printf( "DEBUG: %s\n" , cmd ) ; 
	system( cmd ) ; 
	printf( "running grep... \n" ) ; 
	system( "grep '>' tmp5 > tmp6" ) ; // tmp6 : hit headers  
}

void cleanUp () 
{
	system( "rm ./tmp*" ) ; 
}

