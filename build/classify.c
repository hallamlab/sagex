
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
void classify ( char **sag , int sagN , char **sagNames , char **gm , int gmN , char **gmNames , double alpha , double beta , int threads , double eps , int minLength , int maxIter , int chopSize , int overlap , int proportion, int k , int fixK , double malahanobisMultiple , int verbose , char *kmerFreq , char *kmerPCA, int **out, char *output , int lcsCut ) 
{
	int subDim = 3 ; 
	int cols = 256 ; 
	int rows = gmN + sagN ; 
	int tmp ; 
	int *names = (int*) malloc( gmN * sizeof(int) ) ; // never freed 
	int *sagNamesIdx = (int*) malloc( sagN * sizeof(int) ) ; // never freed 
	double eigenEps = eps * 0.1 ;  
	
	int i , j ; 
	
	/*
	int  *gmSubSet = NULL ; 
	size_t total = 0 ;  
	if( lcsCut > 0 ) 
	{
		if( verbose > 0 ) 
			fprintf( stderr , "Constructing kmer lookup table\n" ) ; 
		
		// get key information  
		int hashKeySize , preKeySize , lastKeySize ; 
		sHashGetKeySize ( lcsCut , &hashKeySize , &preKeySize , &lastKeySize ) ; 
		
		// evaluate SAG entries for entrance into the dictionary 
		size_t *lengths = (size_t*) malloc( sagN * sizeof(size_t) ) ; 
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
		
		gmSubSet = (int*) malloc( gmN * sizeof(int) ) ; 
		
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
		
		// Moved to after Covariance calculation 
		// gmN = total ; // NOTICE UPDATE OF ESSENTIAL VARIABLE 
		// extract subset // Moved to after Covariance calculation  
		for( i = 0 ; i < total ; i++ ) 
		{ 
			gm[i] = gm[ gmSubSet[i] ] ; 
			gmNames[i] = gmNames[ gmSubSet[i] ] ; 
		}  
		
		free( lengths ) ; 
		free( dictionary ) ; 
		free( idx ) ; 
		// free( gmSubSet ) ; 
		free( hashTmp ) ; 
	} 
	*/ 
	
	if( verbose > 0 ) 
		fprintf( stderr , "Calculating kmers for SAG\n" ) ;  
	// Calculate kmer matrix for the SAG 
	int *sagKmers_int = NULL ;  
	// posixCounter( sag , sagN , minLength , threads , chopSize , overlap , &sagKmers_int , &tmp , &sagNamesIdx ) ; 
	int sagKmersN ; 
	countKmers ( sag , sagN , minLength , chopSize , overlap , threads , &sagNamesIdx , &sagKmers_int , &sagKmersN ) ; 
	
//	// reduce sag 
//	for( i = 0 ; i < tmp ; i++ ) // TODO don't do this, you'll destroy valuable information  
//	{ 
//		sag[i] = sag[ sagNamesIdx[i] ] ; 
//		sagNames[i] = sagNames[ sagNamesIdx[i] ] ; 
//	} 
//	sagN = tmp ; // NOTICE CHANGE OF sagN 
	double *sagKmers = (double*) malloc( cols * sagKmersN * sizeof(double) ) ; 
	intToDoubleMat( sagKmers_int , &sagKmersN , &cols , sagKmers ) ;  
	
	if( verbose > 0 )
		fprintf( stderr , "Calculating kmers for Metagenome\n" ) ;
	// Calculate kmer matrix for the gm 
	int *gmKmers_int = NULL ;  
	// posixCounter( gm , gmN , minLength , threads , chopSize , overlap , &gmKmers_int , &tmp , &names ) ;  
	int gmKmersN ; 
	countKmers ( gm , gmN , minLength , chopSize , overlap , threads , &names , &gmKmers_int , &gmKmersN ) ; 
	
//	// reduce gm  
//	for( i = 0 < tmp ; i++ ) // TODO stop, same as above  
//	{ 
//		gm[i] = gm[ names[i] ] ; 
//		gmNames[i] = gmNames[ names[i] ] ; 
//	} 
//	gmN = tmp ; // NOTICE CHANGE OF gmN !!! 
	double *gmKmers = (double*) malloc( cols * gmKmersN * sizeof(double) ) ; 
	intToDoubleMat( gmKmers_int , &gmKmersN , &cols , gmKmers ) ; 
	
	rows = gmKmersN + sagKmersN ; 
	
	// Append matrices 
	double *bigK = (double*) malloc( cols * rows * sizeof(double) ) ; 
	appendRows ( gmKmers , sagKmers , &gmKmersN , &sagKmersN , &cols , bigK ) ; 
	
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
	
	if( kmerFreq != NULL ) // report kmers  
	{
        	FILE *kmerFile;
        	kmerFile = fopen(kmerFreq, "w");
		for( i = 0 ; i < gmKmersN ; i++ )   
		{
			fprintf(kmerFile, "%s" , gmNames[ names[i] ] ) ; 
			for( j = 0 ; j < cols ; j++ ) 
				fprintf(kmerFile, "\t%e" , bigK[ i + rows * j ] ) ; 
			fprintf(kmerFile, "\n" ) ; 
		}
		for( i = 0 ; i < sagKmersN ; i++ ) 
		{
			fprintf(kmerFile, "%s" , sagNames[ sagNamesIdx[i] ] ) ; 
			for( j = 0 ; j < cols ; j++ ) 
				fprintf(kmerFile, "\t%e" , bigK[ gmN + i + rows * j ] ) ; 
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
	double *means = (double*) malloc( cols * sizeof(double) ) ; 
	double *vars = (double*) malloc( cols * sizeof(double) ) ; 
	colStandardize( bigK , &gmKmersN , &rows , &cols , standardK ) ; // TODO Transform w.r.t. gm only! // CHANGE MADE BUT CAUSES FAILURE OF CONV. IN ANDYW.  
	
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
	corrMat ( standardK , &gmKmersN , &rows , &cols , corr ) ; // Calculate PCA basis using GM only! // TODO I expect the error to come from here   
// fprintf( stderr , "DEBUG 3: %f\n" , *corr ) ; 
// double *TMP = (double*) malloc( cols * cols * sizeof(double) ) ; chol( corr , &cols , TMP ) ; double dett = 1.0 ; for( i = 0 ; i < cols ; i++ ){ dett *= TMP[ i + cols * i ] ; } fprintf( stderr , "DEBUG 4: Approximate determinant: %e\n" , dett*dett ) ; free( TMP ) ; 
//int ii ; int nans = 0 ; for( ii = 0 ; ii < cols * cols ; ii++ ){ if( corr[ii] != corr[ii] ) nans++ ; } fprintf( stderr , "DEBUG 4: nans: %i\n" , nans ) ; 
// double dett ; det( corr , &cols , &dett ) ; fprintf( stderr , "DEBUG: det: %e\n" , dett ) ; 
	
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
// double tmpSum = 0.0 ; for( i = 0 ; i < cols ; i++ ) tmpSum += eigVals[i] ;
// fprintf( stderr , "DEBUG 5: sum: %e" , tmpSum ) ; for( i = 0 ; i < cols ; i++ ) fprintf( stderr , " %e" , eigVals[i] ) ; fprintf( stderr , "\n" ) ; 
// fprintf( stderr , "DEBUG 6: %f\n" , *eigVecs ) ; 
	
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
	double *sagC = (double*) malloc( sagKmersN * subDim * sizeof(double) ) ; 
	double *gmC = (double*) malloc( gmKmersN * subDim * sizeof(double) ) ; 
	deappend ( projection , &rows , &subDim , &gmKmersN , gmC , sagC ) ; 
	
	// run identity filter  
	int *idenHits = NULL ; // If indentity filter runs, this will contain the indices of all contigs which pass the identity filter 
	size_t idenHitsN = 0 ;  
	if( lcsCut > 0 ) 
        {  
                idenHits = (int*) malloc( gmN * sizeof(int) ) ; 
                identityFilter ( sag , sagN , gm , gmN , lcsCut , verbose , minLength , idenHits , &idenHitsN ) ; 
        }
	
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
	
	if( k > 1 ) // attempt to fit GMM  
	{ 
		tSAG = (double*) malloc( subDim * sagKmersN * sizeof(double) ) ;
		p = (double*) malloc( k * sizeof(double) ) ; 
		transpose ( sagC , &sagKmersN , &subDim , tSAG ) ; 
		int err = gmmInit ( tSAG , &sagKmersN , &subDim , &k , p , &eps , &threads ,  mu , cov , &maxIter ) ; 
		if( err >= 0 ) 
		{ 
			err = fitMixture ( tSAG , &sagKmersN , &subDim , &k , &eps , p , mu , cov , &maxIter , &threads ) ; 
		} 
		if( err < 0 ) // error! 
		{ 
			if( fixK > 0 ) // cannot change k ! Abort ! 
			{ 
				if( verbose > 0 ) 
					fprintf( stderr , "ERROR: Gaussian Mixture Model failed to fit! Consider reducing -k or deactivating -K\n" ) ; 
				return ; 
			} 
			else 
			{ 
				while( k > 1 && err < 0 ) 
				{ 
					k-- ; 
					if( k > 1 ) 
					{ 
						err = gmmInit ( tSAG , &sagKmersN , &subDim , &k , p , &eps , &threads ,  mu , cov , &maxIter ) ; 
						if( err >= 0 ) 
						{ 
							err = fitMixture ( tSAG , &sagKmersN , &subDim , &k , &eps , p , mu , cov , &maxIter , &threads ) ; 
						} 
					} 
				} 
				if( verbose > 0 ) 
					fprintf( stderr , "WARNING: Argument -k reduced to %i\n" , k ) ; 
			} 
		} 
	} 
	
	if( k == 1 ) 
		covMat ( sagC , &sagKmersN , &subDim , mu , cov ) ; 
	else
	{ 
		// tSAG = (double*) malloc( subDim * sagKmersN * sizeof(double) ) ; 
		// p = (double*) malloc( k * sizeof(double) ) ; 
		// transpose ( sagC , &sagKmersN , &subDim , tSAG ) ; 
		// fitMixture ( tSAG , &sagKmersN , &subDim , &k , &eps , p , mu , cov , &maxIter , &threads ) ; 
		/*
                p[0] = 0.6 ; p[1] = 0.4 ; 
                
                mu[0] = 0.4 ; mu[1] = -0.16 ; mu[2] = 0.315 ; 
                mu[3] = 0.323 ; mu[4] = 0.72 ; mu[5] = -0.0165 ; 
                
                cov[0] = 0.315 ; cov[1] = -0.016 ; cov[2] = -0.24 ; 
                cov[3] = -0.016 ; cov[4] = 0.4266 ; cov[5] = -0.04 ; 
                cov[6] = -0.24 ; cov[7] = -0.04 ; cov[8] = 0.4104 ; 
                
                cov[9] = 0.165 ; cov[10] = -0.018 ; cov[11] = -0.038 ; 
                cov[12] = -0.018 ; cov[13] = 0.027 ; cov[14] = 0.01 ; 
                cov[15] = -0.038 ; cov[16] = 0.01 ; cov[17] = 0.05 ; 
		*/ 
	}
// fprintf( stderr , "DEBUG, mu:\n" ) ; for( i = 0 ; i < subDim ; i++ ){ fprintf( stderr , "%f\n" , mu[i] ) ; } ; fprintf( stderr , "DEBUG, sig:\n" ) ; for( i=0 ; i < subDim ; i++ ){ for( j = 0 ; j < subDim ; j++ ) fprintf( stderr , "%f " , cov[ i + subDim * j ] ) ; fprintf( stderr , "\n" ) ; }
	
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
		fprintf( stderr , "Calculating square root matrix\n" ) ; // TODO move this to the if( k < 2 ) conditional ahead  
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
	matProd ( covEigVecs , diagMat , &subDim , &subDim , &subDim , sqRtMat ) ; // TODO replace with something that doesn't multiply by O(n^2) zeros!   
	
	if( verbose > 0 )
		fprintf( stderr , "Calculating distance statistics\n" ) ;
	// Calculate distance statistics 
	double *diffMat = NULL ; 
	double *statMat = NULL ; 
	int *kmerStatus = (int*) malloc( gmKmersN * sizeof(int) ) ; // 0 : out , 1 : in  
	if( k < 2 ) 
	{ 
		statMat = (double*) malloc( subDim * gmKmersN * sizeof(double) ) ; 
		diffMat = (double*) malloc( subDim * gmKmersN * sizeof(double) ) ; // in R^{ subDim X gmN }  
		
		for( i = 0 ; i < gmKmersN ; i++ ) 
		{
			for( j = 0 ; j < subDim ; j++ ) 
				diffMat[ j + subDim * i ] = (gmC[ i + gmKmersN * j ] - mu[j]) ; 
				// diffMat[ j + subDim * i ] = gmC[ i + gmN * j ] - mu[j] ; 
		} 
		matProd ( sqRtMat , diffMat , &subDim , &subDim , &gmKmersN , statMat ) ;   
		
		// Classify each gm sequence as IN the SAG or not 
		double dtmp ; 
		for( i = 0 ; i < gmKmersN ; i++ ) 
		{
			dtmp = 0.0 ; 
			for( j = 0 ; j < subDim ; j++ ) 
				dtmp += statMat[ j + subDim * i ] * statMat[ j + subDim * i ] ; 
			if( dtmp > malahanobisMultiple*cut )  
				kmerStatus[i] = 0 ; 
			else
				kmerStatus[i] = 1 ;  
		}
	}
	else
	{ 
		int l ; 
		double *invSqrtSig = (double*) malloc( subDim * subDim * k * sizeof(double)  ) ; 
		double *tmpEigMat = (double*) malloc( subDim * subDim * sizeof(double) ) ; 
		double *tmpEigVec = (double*) malloc( subDim * sizeof(double) ) ; 
		double *tmpX = (double*) malloc( subDim * sizeof(double) ) ; 
		for( i = 0 ; i < k ; i++ ) // calc invSqrtSig  
		{ 
			psdEig ( &cov[ subDim * subDim * i ] , &subDim , &eps , tmpEigMat , tmpEigVec ) ; 
			for( j = 0 ; j < subDim ; j++ ) 
			{ 
				for( l = 0 ; l < subDim ; l++ ) 
					invSqrtSig[ j + subDim * l + subDim * subDim * i ] = tmpEigMat[ j + subDim * l ] / sqrt( tmpEigVec[ l ] ) ; 
			} 
		} 
		
		double min , dtmp ; 
		for( i = 0 ; i < gmKmersN ; i++ ) // cycle thru metagenomic points for classification 
		{ 
			/* 
			for( j = 0 ; j < k ; j++ ) // cycle thru gaussians 
			{
				dtmp = 0.0 ; 
				for( l = 0 ; l < subDim ; l++ ) 
					dtmp += pow(fabs(gmC[ i + gmKmersN * l ] - mu[ l + subDim * j ]),2.0) ; 
				dtmp = sqrt(dtmp) ; 
				if( j == 0 ) 
					min = dtmp ; 
				else if ( min > dtmp ) 
					min = dtmp ; 
			} 
			*/ 
			for( j = 0 ; j < subDim ; j++ ) 
				tmpX[j] = gmC[ i + gmKmersN * j ] ; 
			min = calcStat ( tmpX , &k , &subDim , mu , invSqrtSig ) ; 
			
			if( min > malahanobisMultiple*cut ) 
				kmerStatus[i] = 0 ; 
			else
				kmerStatus[i] = 1 ; 
		} 
		free( invSqrtSig ) ; 
		free( tmpEigMat ) ; 
		free( tmpEigVec ) ; 
		free( tmpX ) ; 
	}
	
	if( verbose > 0 )
		fprintf( stderr , "Classifying\n" ) ;
	
	int *categoryTotal = NULL ; 
	int *categoryCount = NULL ; 
	/* 
	if( kmerPCA != NULL ) // only report kmer PCA  
	{
        	FILE *pcaFile;
        	pcaFile = fopen(kmerPCA, "w");
		for( i = 0 ; i < subDim ; i++ ) 
			fprintf(pcaFile, "PC_%i\t" , i ) ; 
		fprintf(pcaFile, "kmerStatus\n" ) ; 
		for( i = 0 ; i < gmKmersN ; i++ ) // cycle through genome entries 
		{
			fprintf(pcaFile, "%s\t" , gmNames[ names[i] ] ) ;  
			for( j = 0 ; j < subDim ; j++ ) 
				fprintf(pcaFile, "%e\t" , gmC[ i + gmKmersN * j ] ) ; 
			fprintf(pcaFile, "%i\n" , kmerStatus[i] ) ; 
		}
		for( i = 0 ; i < sagKmersN ; i++ ) 
		{
			fprintf(pcaFile, "%s\t" , sagNames[ sagNamesIdx[i] ] ) ; 
			for( j = 0 ; j < subDim ; j++ ) 
				fprintf(pcaFile, "%e\t" , sagC[ i + sagKmersN * j ] ) ; 
			fprintf(pcaFile, "2\n" ) ; 
		}
        	fclose(pcaFile); 
	}
	*/
	
	// report hits as fasta  
	categoryTotal = (int*) malloc( gmN * sizeof(int) ) ; 
	categoryCount = (int*) malloc( gmN * sizeof(int) ) ; 
	for( i = 0 ; i < gmN ; i++ )
	{
		categoryTotal[i] = 0 ; 
		categoryCount[i] = 0 ; 
	} 
	for( i = 0 ; i < gmKmersN ; i++ ) 
	{ 
		categoryTotal[ names[i] ] ++ ; 
		if( kmerStatus[i] > 0 ) 
			categoryCount[ names[i] ] ++ ; 
	} 
	
	// reduce hit list via idenHits -- the identity filter 
	int *hits = (int*) malloc( gmN * sizeof(int) ) ; 
	for( i = 0 ; i < gmN ; i++ ) 
	{ 
		hits[i] = 0 ; 
		if( categoryTotal[i] > 0 ) 
		{ 
			if( ((double) categoryCount[i] ) / ((double) categoryTotal[i] ) >= beta ) // minimum acceptable chop percentage to accept a contig is beta-% of its kmers  
			{ 
				// look for i in idenHits 
				if( lcsCut < 1 ) // identity filter is deactivated 
					hits[i] = 1 ; 
				else // identity filter is active   
				{ 
					// TODO use binary search 
					for( j = 0 ; j < idenHitsN && hits[i] == 0 ; j++ ) 
					{ 
						if( idenHits[j] == i ) 
							hits[i] = 1 ;  
					}  
		 		}  
			} 
		} 
	} 
	
        if( kmerPCA != NULL ) // only report kmer PCA  
        {
                FILE *pcaFile;
                pcaFile = fopen(kmerPCA, "w");
                for( i = 0 ; i < subDim ; i++ ) 
                        fprintf(pcaFile, "PC_%i\t" , i ) ; 
                fprintf(pcaFile, "contigStatus\tkmerStatus\n" ) ; 
                for( i = 0 ; i < gmKmersN ; i++ ) // cycle through genome entries 
                {
                        fprintf(pcaFile, "%s\t" , gmNames[ names[i] ] ) ;   
                        for( j = 0 ; j < subDim ; j++ ) 
                                fprintf(pcaFile, "%e\t" , gmC[ i + gmKmersN * j ] ) ; 
                        fprintf(pcaFile, "%i\t%i\n" , hits[ names[i] ] , kmerStatus[i] ) ; 
                }
                for( i = 0 ; i < sagKmersN ; i++ ) 
                {
                        fprintf(pcaFile, "%s\t" , sagNames[ sagNamesIdx[i] ] ) ; 
                        for( j = 0 ; j < subDim ; j++ ) 
                                fprintf(pcaFile, "%e\t" , sagC[ i + sagKmersN * j ] ) ; 
                        fprintf(pcaFile, "2\t2\n" ) ; 
                }
                fclose(pcaFile); 
        }
	
    if (output != NULL)
    {
        FILE *sagexOut;
        sagexOut = fopen(output, "w+");
	for( i = 0 ; i < gmN ; i++ ) 
	{ 
            if( hits[i] > 0 )  
            { 
                fprintf(sagexOut, "%s\n", gmNames[i]);
                fprintf(sagexOut, "%s\n", gm[i]);
            } 
        } 
        fclose(sagexOut);
    }
    else
    {
	for( i = 0 ; i < gmN ; i++ )
	{ 
	    if( hits[i] > 0 ) 
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
	free( hits ) ; 
	free( means ) ; 
	free( vars ) ; 
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
	free( kmerStatus ) ; 
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


