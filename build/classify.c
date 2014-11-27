
#include "classify.h" 

// Classifies contigs as IN a SAG 
// sag : list of standard-length sag sequences 
// sagN : length of the sag array 
// sagNames : a pointer to the list of sag contig names 
// names : integer titles equating only when a metagenomic sequence is from the same contig 
// gm : list of standard-length metagenome / metbag sequences 
// gmNames : a pointer to the list of metagenome contig names 
// alpha : type 1 error rate for classifying an element of gm as IN the SAG 
// beta : proportion of shared names that need to be IN the SAG 
// threads : the desired number of posix threads for this task 
// eps : PCA convergence error parameter 
// minLength : minimum length for contigs  
// maxIter : PCA maximum number of iterations , max iterations in gamma dist integration or max bootstraps in qGMM integration  
// chopsSize : 
// overlap : 
// k : number of gaussians in the gaussian mixture model modelling the SAG 
// verbose : describe process in stderr if > 0 
// kmerFreq : writes kmers to file specified 
// out : a pointer to be allocated with the int-names of contigs which have been classified as IN the SAG 
void classify ( char **sag , int sagN , char **sagNames , char **gm , int gmN , char **gmNames , double alpha , double beta , int threads , double eps , int minLength , int maxIter , int chopSize , int overlap , int proportion, int k , int verbose , char *kmerFreq , char *kmerPCA, int **out, char *output , int lcsCut ) 
{
	int subDim = 3 ; 
	int cols = 256 ; 
	int rows = gmN + sagN ; 
	int tmp ; 
	int *names = NULL ; // never freed 
	int *sagNamesIdx = NULL ; // never freed 
	double eigenEps = eps * 0.0001 ;  
	
	int i , j ; 
	
	 
	if( lcsCut > 0 ) 
	{
		if( verbose > 0 ) 
			fprintf( stderr , "Constructing kmer lookup table\n" ) ; 
		
		// get key information  
		int hashKeySize , preKeySize , lastKeySize ; 
		sHashGetKeySize ( lcsCut , &hashKeySize , &preKeySize , &lastKeySize ) ; 
		
		// evaluate SAG entries for entrance into the dictionary 
		size_t *lengths = (size_t*) malloc( sagN * sizeof(size_t) ) ; 
		size_t total = 0 ; 
		size_t tmpT ; 
		int listN = 0 ; 
		for( i = 0 ; i < sagN ; i++ ) 
		{
			tmpT = strlen( sag[i] ) ; 
			lengths[i] = tmpT ; 
			if( tmpT > minLength ) 
			{
				total += tmpT ; 
				listN += tmpT - lcsCut + 1 ; 
			}
		}
		ull *dictionary = (ull*) malloc( hashKeySize * total * sizeof(ull) ) ;  
		total = 0 ;  
		for( i = 0 ; i < sagN ; i++ ) 
		{ 
			if( lengths[i] > minLength ) 
			{ 
				// record entries 
				sHashes ( sag[i] , hashKeySize , &(dictionary[total]) , preKeySize , lastKeySize , NULL , NULL , 0 ) ; 
				total += (lengths[i] - lcsCut + 1) * hashKeySize ;  
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
		
		int *gmSubSet = (int*) malloc( gmN * sizeof(int) ) ; 
		
		free( lengths ) ; 
		lengths = (size_t*) malloc( gmN * sizeof(size_t) ) ; 
		tmp = 0 ; 
		for( i = 0 ; i < gmN ; i++ ) 
		{
			lengths[i] = strlen( gm[i] ) ; 
			if( lengths[i] > tmp ) 
				tmp = lengths[i] ; 
		} 
		ull *hashTmp = (ull*) malloc( hashKeySize * tmp * sizeof(ull) ) ; 
		
		if( verbose > 0 ) 
			fprintf( stderr , "Performing kmer lookups\n" ) ; 
		
		total = 0 ; 
		for( i = 0 ; i < gmN ; i++ ) 
		{ 
			if( lengths[i] > minLength ) 
			{ 
				tmp = sHashes ( gm[i] , hashKeySize , hashTmp , preKeySize , lastKeySize , dictionary , idx , listN ) ; 
				if( tmp > 0 ) 
				{ 
					gmSubSet[total] = i ; 
					total++ ; 
				} 
			} 
			if( verbose > 0 ) 
				fprintf( stderr , "\rGM: %f%%" , 100.0f*((float) i)/((float) gmN) ) ; 
		} 
		if( verbose > 0 ) 
			fprintf( stderr , "\rGM: 100.0%%          \nExtracted subset of size %lu\n" , total ) ; 
		gmN = total ; // NOTICE UPDATE OF ESSENTIAL VARIABLE 
		
		// extract subset 
		for( i = 0 ; i < total ; i++ ) 
		{ 
			gm[i] = gm[ gmSubSet[i] ] ; 
			gmNames[i] = gmNames[ gmSubSet[i] ] ; 
		} 
		
		free( lengths ) ; 
		free( dictionary ) ; 
		free( idx ) ; 
		free( gmSubSet ) ; 
		free( hashTmp ) ; 
	}
	
	if( verbose > 0 ) 
		fprintf( stderr , "Calculating kmers for SAG\n" ) ;  
	// Calculate kmer matrix for the SAG 
	int *sagKmers_int = NULL ; // = (int*) malloc( cols * sagN * sizeof(int) ) ; 
	// posixCounter( sag , sagN , minLength , threads , chopSize , overlap , &sagKmers_int , &tmp , &sagNamesIdx ) ; 
	countKmers ( sag , sagN , minLength , chopSize , overlap , threads , &sagNamesIdx , &sagKmers_int , &tmp ) ; 
	sagN = tmp ; // NOTICE CHANGE OF sagN 
	double *sagKmers = (double*) malloc( cols * sagN * sizeof(double) ) ; 
	intToDoubleMat( sagKmers_int , &sagN , &cols , sagKmers ) ; 
	
	if( verbose > 0 )
		fprintf( stderr , "Calculating kmers for Metagenome\n" ) ;
	// Calculate kmer matrix for the gm 
	int *gmKmers_int = NULL ; // (int*) malloc( cols * gmN * sizeof(int) ) ; 
	// posixCounter( gm , gmN , minLength , threads , chopSize , overlap , &gmKmers_int , &tmp , &names ) ;  
	countKmers ( gm , gmN , minLength , chopSize , overlap , threads , &names , &gmKmers_int , &tmp ) ;
	gmN = tmp ; // NOTICE CHANGE OF  gmN  
	double *gmKmers = (double*) malloc( cols * gmN * sizeof(double) ) ; 
	
	rows = gmN + sagN ;  
	
	/*
	int k ; 
	for( k = 0 ; k < gmN ; k++ ) 
		printf( "%s\n" , gm[k] ) ; 
	*/
	
	/* 
	for( i = 0 ; i < sagN ; i++ ) 
	{ 
		for( j = 0 ; j < cols ; j++ ) 
		{
			// printf( "DEBUG seg, k: %i, l: %i\n" , k , l ) ; 
			fprintf( stderr , "%i " , sagKmers_int[i + sagN * j] ) ;  
		}
		fprintf( stderr , "\n" ) ; 
	}*/ 
	
	intToDoubleMat( gmKmers_int , &gmN , &cols , gmKmers ) ; 
	
	// Append matrices 
	double *bigK = (double*) malloc( cols * rows * sizeof(double) ) ; 
	appendRows ( gmKmers , sagKmers , &gmN , &sagN , &cols , bigK ) ; 
	
	// convert to proportions, if requested 
        if( proportion == 1 ) 
        {   
                double sum ; 
                for( i = 0 ; i < rows ; i++ ) 
                {   
                        sum = 0.0 ; 
                        for( j = 0 ; j < cols ; j++ ) 
                                sum += bigK[ i + rows * j ] ; 
                        for( j = 0 ; j < cols ; j++ ) 
                                bigK[ i + rows * j ] = bigK[ i + rows * j ] / sum ; 
                }   
        }
	
	if( kmerFreq != NULL ) // report kmers and quit 
	{
        FILE *kmerFile;
        kmerFile = fopen(kmerFreq, "w");
		for( i = 0 ; i < gmN ; i++ )   
		{
			fprintf(kmerFile, "%s" , gmNames[ names[i] ] ) ; 
			for( j = 0 ; j < cols ; j++ ) 
				fprintf(kmerFile, "\t%e" , bigK[ i + (gmN + sagN) * j ] ) ; 
			fprintf(kmerFile, "\n" ) ; 
		}
		for( i = 0 ; i < sagN ; i++ ) 
		{
			fprintf(kmerFile, "%s" , sagNames[ sagNamesIdx[i] ] ) ; 
			for( j = 0 ; j < cols ; j++ ) 
				fprintf(kmerFile, "\t%e" , bigK[ gmN + i + (gmN + sagN) * j ] ) ; 
			fprintf(kmerFile, "\n" ) ; 
		}
        fclose(kmerFile);
	}
	
	/*
	for( i = 0 ; i < rows ; i++ ) 
	{
		for( j = 0 ; j < cols ; j++ ) 
			fprintf( stderr , "%f " , bigK[ i + rows * j ] ) ; 
		fprintf( stderr , "\n" ) ; 
	}*/
	
	// standardizes dimensions 
	double *standardK = (double*) malloc( cols * rows * sizeof(double) ) ; 
	colStandardize( bigK , &rows , &cols , standardK ) ; 
	
	/*
	{int k , l ; 
	for( k = 0 ; k < rows ; k++ ) 
	{
		for( l = 0 ; l < cols ; l++ ) 
			fprintf( stderr , "%f " , standardK[ k + rows * l ] ) ; 
		fprintf( stderr , "\n" ) ; 
	}}*/
	
	if( verbose > 0 ) 
		fprintf( stderr , "Calculating correlation matrix\n" ) ;
	// calculate correlation matrix 
	double *corr = (double*) malloc( cols * cols * sizeof(double) ) ; 
	corrMat ( standardK , &rows , &cols , corr ) ; 
	
	/*
	for( i = 0 ; i < cols ; i++ ) 
	{
		for( j = 0 ; j < cols ; j++ ) 
			fprintf( stderr , "%f " , corr[i + cols * j] ) ; 
		fprintf( stderr , "\n" ) ; 
	}
	*/
	
	if( verbose > 0 ) 
		fprintf( stderr , "Calculating eigen vectors\n" ) ; 
	// Calculate eigen space 
	double *eigVecs = (double*) malloc( cols * cols * sizeof(double) ) ; 
	double *eigVals = (double*) malloc( cols * sizeof(double) ) ; 	
	// powerIteration ( corr , &cols , &subDim , &eigenEps , eigVals , eigVecs , &threads ) ;  
	psdEig ( corr , &cols , &eigenEps , eigVecs , eigVals ) ; 
	
	/*
	printf( "Eigen values: " ) ; 
	for( i = 0 ; i < subDim ; i++ ) 
		fprintf( stderr , "%f " , eigVals[i] ) ; 
	fprintf( stderr , "\n" ) ; 
	*/
	/*
	for( i = 0 ; i < cols ; i++ ) 
	{
		for( j = 0 ; j < cols ; j++ ) 
			fprintf( stderr , "%f " , eigVecs[ i + cols * j ] ) ; 
		fprintf( stderr , "\n" ) ; 
	}*/
	
	if( verbose > 0 ) 
		fprintf( stderr , "Projecting into principle component subspace\n" ) ;
	// Project into first three principal components 
	double *projection = (double*) malloc( rows * subDim * sizeof(double) ) ; 
	matProd ( standardK , eigVecs , &rows , &cols , &subDim , projection ) ; 
	
	// Extract sag and gm clusters 
	double *sagC = (double*) malloc( sagN * subDim * sizeof(double) ) ; 
	double *gmC = (double*) malloc( gmN * subDim * sizeof(double) ) ; 
	deappend ( projection , &rows , &subDim , &gmN , gmC , sagC ) ; 
	
	/*
	int ii , jj ; 
	for( ii = 0 ; ii < gmN ; ii++ ) 
	{
		for( jj = 0 ; jj < subDim ; jj++ ) 
			fprintf( stdout , "%f " , gmC[ ii + gmN*jj ] ) ; 
		fprintf( stdout , "\n" ) ; 
	}
	for( ii = 0 ; ii < rows - gmN ; ii++ ) 
	{
		for( jj = 0 ; jj < subDim ; jj++ ) 
			fprintf( stderr , "%f " , sagC[ ii + (rows-gmN)*jj ] ) ; 
		fprintf( stderr , "\n" ) ; 
	}
	*/
	
	if( verbose > 0 ) 
		fprintf( stderr , "Training classifier\n" ) ; 
	// Train the classifier on the SAG 
	double *mu = (double*) malloc( subDim * k * sizeof(double) ) ; 
	double *cov = (double*) malloc( subDim * subDim * k * sizeof(double) ) ; 
	double *tSAG = NULL ; 
	double *p = NULL ; 
	if( k == 1 ) 
		covMat ( sagC , &sagN , &subDim , mu , cov ) ; 
	else
	{ 
		tSAG = (double*) malloc( subDim * sagN * sizeof(double) ) ; 
		p = (double*) malloc( k * sizeof(double) ) ; 
		transpose ( sagC , &sagN , &subDim , tSAG ) ; 
		fitMixture ( tSAG , &sagN , &subDim , &k , &eps , p , mu , cov , &maxIter , &threads ) ; 
	}
	
	if( verbose > 0 )
		fprintf( stderr , "Designing cut-off\n" ) ;
	// Calc cut-off w.r.t. ref. dist. 
	double cut ; 
	if( k < 2 ) 
	{
		double cumulativeP = 1.0 - alpha ; 
		double rate = 2.0 ; 
		double shape = ((double) subDim ) / 2.0 ;  
		dyadicInvGam ( &shape , &rate , &cumulativeP , &eps , &maxIter , &cut ) ;  
	}
	else
	{   
// int i , j ; fprintf( stderr , "p: " ) ; for( i = 0 ; i < k ; i++ ) { fprintf( stderr , " %e" , p[i] ) ; } fprintf( stderr , "\n" ) ; 
// fprintf( stderr , "mu:\n" ) ; for( i = 0 ; i < subDim ; i++ ){ for( j = 0 ; j < k ; j++ ) fprintf( stderr , "%e " , mu[ i + subDim * j ] ) ; fprintf( stderr , "\n" ) ; } 
// int ll; for( ll = 0 ; ll < k ; ll++ ){ fprintf( stderr , "Cov %i:\n" , ll ) ; for( i = 0 ; i < subDim ; i++ ){ for( j = 0 ; j < subDim ; j++ ) fprintf( stderr , "%e " , cov[ i + subDim * j + subDim * subDim * ll ] ) ; fprintf( stderr , "\n" ) ; } }
		double quantile = 1.0 - alpha ; 
		getGMMQuantile ( &quantile , p , mu , cov , &maxIter , &subDim , &k , &eps , &threads , &cut ) ; 
	}
// fprintf( stderr , "Cut: %e\n" , cut ) ; 
	
	if( verbose > 0 )
		fprintf( stderr , "Calculating square root matrix\n" ) ;
	// Calculate square root matrix of covariance 
	double *sqRtMat = (double*) malloc( subDim * subDim * sizeof(double) ) ; 
	double *covEigVecs = (double*) malloc( subDim * subDim * sizeof(double) ) ; 
	double *covEigVals = (double*) malloc( subDim * sizeof(double) ) ; 
	// qrEig ( cov , &subDim , &eps, &minIter , &maxIter , covEigVals , covEigVecs ) ; 
	psdEig( cov , &subDim , &eps , covEigVecs , covEigVals ) ; 
	double *diagMat = (double*) malloc( subDim * subDim * sizeof(double) ) ; 
	for( i = 0 ; i < subDim ; i++ )
	{
		for( j = 0 ; j < subDim ; j++ ) 
		{
			if( i == j ) 
				diagMat[i + subDim * j] = 1.0 / sqrt( covEigVals[i] ) ; // Actually inverse of the square root matrix 
			else
				diagMat[i + subDim * j] = 0.0 ; 
		}
	}
	matProd ( covEigVecs , diagMat , &subDim , &subDim , &subDim , sqRtMat ) ; 
	
	if( verbose > 0 )
		fprintf( stderr , "Calculating distance statistics\n" ) ;
	// Calculate distance statistics 
	double *diffMat = NULL ; 
	double *statMat = NULL ; 
	int *status = (int*) malloc( gmN * sizeof(int) ) ; // 0 : out , 1 : in 
	if( k < 2 ) 
	{
		diffMat = (double*) malloc( subDim * gmN * sizeof(double) ) ; // in R^{ subDim X gmN }  
		statMat = (double*) malloc( subDim * gmN * sizeof(double) ) ; 
		for( i = 0 ; i < gmN ; i++ ) 
		{
			for( j = 0 ; j < subDim ; j++ ) 
				diffMat[ j + subDim * i ] = (gmC[ i + gmN * j ] - mu[j]) ; 
				// diffMat[ j + subDim * i ] = gmC[ i + gmN * j ] - mu[j] ; 
		}
		matProd ( sqRtMat , diffMat , &subDim , &subDim , &gmN , statMat ) ;   
		
		// Classify each gm sequence as IN the SAG or not 
		double dtmp ; 
		for( i = 0 ; i < gmN ; i++ ) 
		{
			dtmp = 0.0 ; 
			for( j = 0 ; j < subDim ; j++ ) 
				dtmp += statMat[ j + subDim * i ] * statMat[ j + subDim * i ] ; 
			if( dtmp > cut ) 
				status[i] = 0 ; 
			else
				status[i] = 1 ;  
		}
	}
	else
	{
		int l ; 
		double min , dtmp ; 
		for( i = 0 ; i < gmN ; i++ ) // cycle thru metagenomic points for classification 
		{
			for( j = 0 ; j < k ; j++ ) // cycle thru gaussians 
			{
				dtmp = 0.0 ; 
				for( l = 0 ; l < subDim ; l++ ) 
					dtmp += pow(fabs(gmC[ i + gmN * l ] - mu[ l + subDim * j ]),2.0) ; 
				dtmp = sqrt(dtmp) ; 
				if( j == 0 ) 
					min = dtmp ; 
				else if ( min > dtmp ) 
					min = dtmp ; 
			}
			if( min > cut ) 
				status[i] = 0 ; 
			else
				status[i] = 1 ; 
		}
	}
	
	if( verbose > 0 )
		fprintf( stderr , "Classifying\n" ) ;
	// Count proportion of chops in the SAG and compare to beta before storing in out 
	int max = 0 ; 
	for( i = 0 ; i < gmN ; i++ ) 
	{
		if( names[i] > max ) 
			max = names[i] ; 
	}
	
	int *categoryTotal = NULL ; 
	int *categoryCount = NULL ; 
	if( kmerPCA != NULL ) // only report kmer PCA  
	{
        FILE *pcaFile;
        pcaFile = fopen(kmerPCA, "w");
		for( i = 0 ; i < subDim ; i++ ) 
			fprintf(pcaFile, "PC_%i\t" , i ) ; 
		fprintf(pcaFile, "status\n" ) ; 
		for( i = 0 ; i < gmN ; i++ ) // cycle through genome entries 
		{
			fprintf(pcaFile, "%s\t" , gmNames[ names[i] ] ) ;  
			for( j = 0 ; j < subDim ; j++ ) 
				fprintf(pcaFile, "%e\t" , gmC[ i + gmN * j ] ) ; 
			fprintf(pcaFile, "%i\n" , status[i] ) ; 
		}
		for( i = 0 ; i < sagN ; i++ ) 
		{
			fprintf(pcaFile, "%s\t" , sagNames[ sagNamesIdx[i] ] ) ; 
			for( j = 0 ; j < subDim ; j++ ) 
				fprintf(pcaFile, "%e\t" , sagC[ i + sagN * j ] ) ; 
			fprintf(pcaFile, "2\n" ) ; 
		}
        fclose(pcaFile);
	}
// report hits as fasta  
//		*out = (int*) malloc( max * sizeof(int) ) ; 
	categoryTotal = (int*) malloc( max * sizeof(int) ) ; 
	categoryCount = (int*) malloc( max * sizeof(int) ) ; 
	for( i = 0 ; i < max ; i++ )
	{
		categoryTotal[i] = 0 ; 
		categoryCount[i] = 0 ; 
	}
	for( i = 0 ; i < gmN ; i++ ) 
	{
		categoryTotal[ names[i] ] ++ ; 
		if( status[i] > 0 ) 
			categoryCount[ names[i] ] ++ ; 
	}
    if (output != NULL)
    {
        FILE *sagexOut;
        sagexOut = fopen(output, "w+");
	    for( i = 0 ; i < max ; i++ ) 
	    {
            if( ((double) categoryCount[i] ) / ((double) categoryTotal[i] ) >= beta )// allow for 100% 
            {
                fprintf(sagexOut, "%s\n", gmNames[i]);
                fprintf(sagexOut, "%s\n", gm[i]);
            }
        }
        fclose(sagexOut);
    }
    else
    {
	    for( i = 0 ; i < max ; i++ )
	    {
            if( ((double) categoryCount[i] ) / ((double) categoryTotal[i] ) >= beta )// allow for 100% 
            {
                fprintf(stdout, "%s\n", gmNames[i]);
                fprintf(stdout, "%s\n", gm[i]);
            }
        }
    }
                
/*
			if( categoryTotal[i] == 0 ) 
			{
				(*out)[i] = 0 ; // if not found, it is definitely not in the SAG 
			}
			else
			{
				if( ((double) categoryCount[i] ) / ((double) categoryTotal[i] ) >= beta ) // allow for 100% 
					(*out)[i] = 1 ; // in the SAG 
				else
					(*out)[i] = 0 ; // not in the SAG 
            }
*/

	// free unused memory  
	free( sagKmers_int ) ; 
	free( gmKmers_int ) ; 
	free( sagKmers ) ; 
	free( gmKmers ) ; 
	free( bigK ) ; 
	free( standardK ) ; 
	free( corr ) ; 
	free( eigVals ) ; 
	free( eigVecs ) ; 
	free( projection ) ; 
	free( gmC ) ; 
	free( sagC ) ; 
	free( mu ) ; 
	free( cov ) ; 
	free( sqRtMat ) ; 
	free( covEigVecs ) ; 
	free( covEigVals ) ; 
	free( diagMat ) ; 
	free( status ) ; 
	if( categoryTotal != NULL ) 
		free( categoryTotal ) ; 
	if( categoryCount != NULL ) 
		free( categoryCount ) ; 
	if( tSAG != NULL ) 
		free( tSAG ) ; 
	if( p != NULL ) 
		free( p ) ; 
	if( diffMat != NULL ) 
		free( diffMat ) ; 
	if( statMat != NULL ) 
		free( statMat ) ; 
}


