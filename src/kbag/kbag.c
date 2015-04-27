
#include "kbag.h" 

typedef unsigned long long ull ; 

int min( int a , int b ) 
{
	return (a > b)? b : a ; 
}

// string hashes are arrays of ull, constant length n 
// returns +1 if a > b, -1 if a < b , 0 if a = b  
int sHashCompare ( ull *a , ull *b , int n ) 
{ 
	int i ; 
	for( i = 0 ; i < n ; i++ ) 
	{
		if( a[i] < b[i] ) 
			return 1 ; 
		if( a[i] > b[i] ) 
			return -1 ; 
	}
	return 0 ; 
} 

// returns index if found, -1 otherwise 
// key : the value to look up 
// list : the memory location of the dictionary 
// idx : the sorted indices of list entries 
// listN : the number of list entries , assumed > 0  
// n : the number of ull's required to store a hash  
int sHashBinarySearch ( ull *key , ull *list , int *idx , int listN , int n ) 
{
	int upper = 0 ; 
	int lower = listN-1 ; 
	int mid = (upper + lower)/2 ; 
	int c ; 
	
	while( upper <= lower ) 
	{
		c = sHashCompare( key , &(list[idx[mid]]) , n ) ; 
		if( c < 0 ) 
		{
			upper = mid + 1 ; 
		}
		else if( c == 0 ) 
		{
			return idx[mid] ; 
		}
		else
		{
			lower = mid - 1 ; 
		}
		mid = (upper + lower)/2 ; 
	} 
	return -1 ; 
}

// TODO this currently allows for duplicate storage. This could result in a RAM explosion! 
// quick sorts a list of sHashes 
// list : storage space for the ulls 
// idx : storage place for the ull IDs 
// listN : number of sHashes in list 
// n : number of ull's per sHash 
void sHashQuickSort ( ull *list , int *idx , int n , int listN ) 
{
	if( listN <= 1 ) 
		return ; 
	int tmp , c ; 
	if( listN == 2 ) 
	{ 
		c = sHashCompare( &(list[idx[0]]) , &(list[idx[1]]) , n ) ; 
		if( c < 0 ) // idx_0 > idx_1 
		{
			tmp = idx[1] ; 
			idx[1] = idx[0] ; 
			idx[0] = tmp ; 
		}
		return ; 
	} 
	int pivot = listN - 1 ; 
	int i ; 
	int j = 0 ; 
	for( i = 0 ; i < listN-1 ; i++ ) 
	{
		c = sHashCompare( &(list[idx[i]]) , &(list[idx[pivot]]) , n ) ; 
		if( c > 0 ) // idx_i > idx_pivot 
		{
			tmp = idx[i] ; 
			idx[i] = idx[j] ; 
			idx[j] = tmp ; 
			j++ ; 			
		}
	} 
	tmp = idx[j] ; 
	idx[j] = idx[pivot] ; 
	idx[pivot] = tmp ;
// printf( "DEBUG j: %i, " , j ) ;  
// for( i = 0 ; i < listN ; i++ ) printf( "%llu " , list[idx[i]] ) ; printf( "\n" ) ; 
	
	sHashQuickSort( list , idx , n , j ) ; 
	sHashQuickSort( list , &(idx[j+1]) , n , listN - j - 1 ) ; 
}

// get number of ull's required to store a key for a kmer of size k  
// k : kmer length 
// n : output , number of ull's required to store a key for a kmer of size k 
// preN : number of characters to represent in the first n-1 hash chunks 
// lastN : number of character to represent in the n-th hash chunk    
void sHashGetKeySize ( int k , int *n , int *preN , int *lastN ) 
{ 
	// calculate the number of characters to be stored per  
	ull space = ULLONG_MAX ; 
	*preN = 0 ; 
	while( space > 3 && *preN < k ) // counting protects against roundoff errors 
	{ 
		space /= 4 ; 
		(*preN)++ ; 
	} 
	
	*n = 1 ; 
	if( *preN == k )  
		*lastN = *preN ; 
	else 
	{
		*n = k / *preN ; 
		if( *preN * *n < k ) 
		{
			*lastN = k - *n * *preN ; 
			(*n)++ ; 
		}
		else
		{
			*lastN = *preN ; 
		}
	}
} 

// defunkt 
// calculate the hash for a string 
// not designed to calculate hashes for all substrings 
// s : the string 
// m : number of characters to store per hash chunk 
// h : where to store the hash 
void sHash ( char *s , int n , int m , ull *h ) 
{
	int i = n-1 ; // h index, notice the big endian  
	int j = 0 ; // characters since last hash chunk 
	
	h[i] = 0 ; 
	ull idx = 0 ; 
	while( *s != '\0' ) 
	{
		switch(*s) 
		{
			case 'A': 
			case 'a': 
				idx = 4*idx + 0 ; 
				break ; 
			case 'T': 
			case 't': 
				idx = 4*idx + 1 ; 
			case 'C': 
			case 'c': 
				idx = 4*idx + 2 ; 
			case 'G': 
			case 'g': 
				idx = 4*idx + 3 ; 
			case 'N': 
			case 'n': 
				break ; // skip unknowns 
			default: 
				fprintf( stderr , "ERROR: sHash, bad character: %c\n" , *s ) ; 
				return ; 
		}
		j++ ; 
		if( j == m ) 
		{
			j = 0 ; 
			h[i] = idx ; 
			idx = 0 ;  
			i-- ; 
		}
		s++ ; 
	}
	if( j < m ) 
		h[i] = idx ; 
}

// calculates the hashes for all substrings of length cut 
// a string of length N will have N - cut hashes 
// s : the string to be have substrings hashed 
// n : the number of ull's required to store a hash of size k  
// dict : if not null, sHashes will perform lookups.  
// sHashes will return 1 if a kmer from s is in dict , 0 otherwise. 
// idx : must not be null if dict is not null. Stores the sorted indices of dict 
// idxN : length of idx 
// sHashes returns -1 on errors  
int sHashes ( char *s , int n , ull *hash , int preN , int lastN , ull *dict , int *idx , int idxN ) 
{ 
	int tmp ; 
	ull preC = 1 ; // pre-chunk size 
	ull lastC = 1 ; // post-chunk size 
	int i ; 
	for( i = 0 ; i < preN ; i++ ) 
		preC *= 4 ; 
	for( i = 0 ; i < lastN ; i++ ) 
		lastC *= 4 ; 
	 
	ull *chunks = (ull*) malloc( n * sizeof(ull) ) ; 
	char **pos = (char**) malloc( n * sizeof(char*) ) ; 
	for( i = 0 ; i < n-1 ; i++ ) 
	{ 
		chunks[i] = 0 ; 
		pos[i] = s + i * preN ; // starting position is a multiple of preN 
	}
	if( n > 1 )  
		pos[n-1] = s + (n-2) * preN + lastN ; 
	else 
		pos[n-1] = s ; 
	
	int j = 0 ; 
	while( *(pos[n-1]) != '\0' ) 
	{ 
		for( i = 0 ; i < n ; i++ ) 
		{ 
			switch( *(pos[i]) ) 
			{ 
				case 'A': 
				case 'a': 
					chunks[i] = 4*chunks[i] + 0 ; 
					break ; 
				case 'T': 
				case 't': 
					chunks[i] = 4*chunks[i] + 1 ; 
					break ; 
				case 'C': 
				case 'c': 
					chunks[i] = 4*chunks[i] + 2 ; 
					break ; 
				case 'G': 
				case 'g': 
					chunks[i] = 4*chunks[i] + 3 ; 
					break ; 
				case 'N': 
				case 'n': 
					break ; // ignore unkowns 
				default: 
					fprintf( stderr , "ERROR: sHashes: invalid character: %c\n" , *(pos[i]) ) ; 
					return -1 ; 
			} 
			(pos[i])++ ; 
		}
		j++ ; 
		if( j >= preN ) 
		{ 
			for( i = 0 ; i < n-1 ; i++ ) 
			{ 
				chunks[i] = chunks[i] % preC ; 
				hash[ i + n*(j-preN) ] = chunks[i] ; 
			} 
			chunks[n-1] = chunks[n-1] % lastC ; 
			hash[ n-1 + n*(j-preN)] = chunks[n-1] ; 
			if( dict != NULL ) 
			{ 
				tmp = sHashBinarySearch ( &(hash[n*(j-preN)]) , dict , idx , idxN , n ) ; 
				if( tmp > -1 ) 
					return 1 ; 
			} 
		}
	}
	
	if( j < preN ) 
		fprintf( stderr , "WARNING: sHashes, contig too short. j: %i, preN: %i\n" , j , preN ) ; 
	
	free( chunks ) ; 
	free( pos ) ; 
	
	return 0 ; 
} 

// Identifies which contigs from gm share at least cut continuous bases with at least one contig from sag 
// sag : the sequences from the sag, length sagN 
// gm : the sequences from the metagenome, length gmN 
// cut : the cutoff, minimum shared bases 
// verbose : > 0 if progress should be printed to stderr 
// minLength : minimum contig length 
// gmSubSet : output, pre-allocated length of gmN  
// gmSubSetN : output, length of gmSubSet used 
void identityFilter ( char **sag , int sagN , char **gm , int gmN , int cut , int verbose , int minLength , int *gmSubSet , size_t *gmSubSetN ) 
{ 
                if( verbose > 0 ) 
                        fprintf( stderr , "Constructing kmer lookup table\n" ) ; 
    
                // get key information  
                int hashKeySize , preKeySize , lastKeySize ; 
                sHashGetKeySize ( cut , &hashKeySize , &preKeySize , &lastKeySize ) ; 
    
                // evaluate SAG entries for entrance into the dictionary 
                size_t *lengths = (size_t*) malloc( sagN * sizeof(size_t) ) ; 
                size_t tmpT ; 
                int listN = 0 ; 
		*gmSubSetN = 0 ; 
		int i ; 
                for( i = 0 ; i < sagN ; i++ ) 
                {   
			if( sag[i] != NULL ) // TODO strings should NEVER be null !!! 
			{ 
                        	tmpT = strlen( sag[i] ) ; 
                        	lengths[i] = tmpT ; 
                        	if( tmpT > (size_t) minLength ) 
                        	{   
                        	        *gmSubSetN += tmpT ; 
                        	        listN += tmpT - cut + 1 ; 
                        	} 
			}   
			else 
				lengths[i] = 0 ; 
                }   
                ull *dictionary = (ull*) malloc( hashKeySize * (*gmSubSetN) * sizeof(ull) ) ;   
                *gmSubSetN = 0 ;   
                for( i = 0 ; i < sagN ; i++ ) 
                {   
                        if( lengths[i] > (size_t) minLength ) 
                        {   
                                // record entries 
                                sHashes ( sag[i] , hashKeySize , &(dictionary[*gmSubSetN]) , preKeySize , lastKeySize , NULL , NULL , 0 ) ; 
                                *gmSubSetN += (lengths[i] - cut + 1) * hashKeySize ;   
                        }
                        if( verbose > 0 )
                                fprintf( stderr , "\rSAG: %f%%" , 100.0f*((float) i)/((float) sagN) ) ;
                }
                if( verbose > 0 )
                        fprintf( stderr , "\rSAG: 100.0%%        \n" ) ;
                int *idx = (int*) malloc( listN * sizeof(int) ) ;
                for( i = 0 ; i < listN ; i++ )
                        idx[i] = hashKeySize * i ; 
		
		                sHashQuickSort ( dictionary , idx , hashKeySize , listN ) ;

                // gmSubSet = (int*) malloc( gmN * sizeof(int) ) ; // now occurs outside of this function 

                free( lengths ) ;
                lengths = (size_t*) malloc( gmN * sizeof(size_t) ) ;
                int tmp = 0 ;
                for( i = 0 ; i < gmN ; i++ )
                { 
			if( gm[i] != NULL ) // Jerry rigged robustness TODO strings should NEVER be null !!! 
			{ 
                        	lengths[i] = strlen( gm[i] ) ;
                        	if( lengths[i] > (size_t) tmp )
                        	        tmp = lengths[i] ; 
			} 
			else 
				lengths[i] = 0 ; 
                } 
                ull *hashTmp = (ull*) malloc( hashKeySize * tmp * sizeof(ull) ) ;

                if( verbose > 0 )
                        fprintf( stderr , "Performing kmer lookups\n" ) ;

                *gmSubSetN = 0 ;
                for( i = 0 ; i < gmN ; i++ )
                {
                        if( lengths[i] > (size_t) minLength )
                        {
                                tmp = sHashes ( gm[i] , hashKeySize , hashTmp , preKeySize , lastKeySize , dictionary , idx , listN ) ;
                                if( tmp > 0 )
                                {
                                        gmSubSet[*gmSubSetN] = i ; 
                                        *gmSubSetN += 1 ;
                                }
                        }
                        if( verbose > 0 )
                                fprintf( stderr , "\rGM: %f%%" , 100.0f*((float) i)/((float) gmN) ) ;
                }
                if( verbose > 0 )
                        fprintf( stderr , "\rGM: 100.0%%          \nExtracted subset of size %lu\n" , *gmSubSetN ) ;	
		
		free( lengths ) ;
                free( dictionary ) ;
                free( idx ) ;
                free( hashTmp ) ;
} 









