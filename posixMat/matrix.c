
#include "matrix.h"

void quickSort ( double *x , int *n , int *idx ) 
{
        if( *n <= 1 ) 
                return ; 
        double tmp ; 
	int tmpIdx ; 
        if( *n == 2 ) 
        {   
                if( x[0] > x[1] ) 
                {   
                        tmp = x[1] ; 
                        x[1] = x[0] ; 
                        x[0] = tmp ; 
			
			if( idx != NULL ) 
			{
				tmpIdx = idx[1] ; 
				idx[1] = idx[0] ; 
				idx[0] = tmpIdx ; 
			}
                }   
                return ; 
        }   
    
        int pivot = *n - 1 ;   
    
        int i ; 
        int j = 0 ; 
        for( i = 0 ; i < *n-1 ; i++ ) 
        {   
                if( x[i] < x[pivot] ) 
                {   
                        tmp = x[i] ; 
                        x[i] = x[j] ; 
                        x[j] = tmp ; 
			
			if( idx != NULL ) 
			{
				tmpIdx = idx[i] ; 
				idx[i] = idx[j] ; 
				idx[j] = tmpIdx ; 
			}
			
			j++ ; 
                }   
        }   
        tmp = x[pivot] ; 
        x[pivot] = x[j] ; 
        x[j] = tmp ; 
	
	if( idx != NULL ) 
	{
		tmpIdx = idx[pivot] ; 
		idx[pivot] = idx[j] ; 
		idx[j] = tmpIdx ; 
	}
	
	pivot = j ; 
    
        i = pivot  ;   
        quickSort( x , &i , idx ) ; 
        i = *n - pivot - 1 ; 
	if( idx != NULL ) 
        	quickSort( &x[pivot+1] , &i , &idx[pivot+1] ) ; 
	else 
		quickSort( &x[pivot+1] , &i , NULL ) ; 
}

// transpose a matrix 
// x : the matrix, also the output space 
// m : number of rows 
// n : number of columns 
// out : space for n X m doubles 
void transpose ( double *x , int *m , int *n , double *out ) 
{
	int i , j ; 
	for( i = 0 ; i < *m ; i++ ) 
	{
		for( j = 0 ; j < *n ; j++ ) 
			out[ j + *n * i ] = x[ i + *m * j ] ; 
	}
}

// subtract matrix b from a 
// n : rows 
// m : cols 
// out : space for n X m doubles 
void subtract( double *a , double *b , int *n , int *m , double *out ) 
{
	int i , j ; 
	for( i = 0 ; i < *n ; i++ ) 
	{
		for( j = 0 ; j < *m ; j++ ) 
			out[ i + *n * j ] = a[ i + *n * j ] - b[ i + *n * j ] ; 
	}
}

// Break a matrix into two matrices by dividing at a row 
// in : a column-major order matrix of n X m entries 
// div : number of rows in out1, from the upper part of 'in'  
// out1 , out2 : pre-allocated space 
void deappend ( double *in , int *n , int *m , int *div , double *out1 , double *out2 ) 
{
	int i , j ; 
	for( i = 0 ; i < *div ; i++ ) 
	{
		for( j = 0 ; j < *m ; j++ ) 
			out1[ i + *div * j ] = in[ i + *n * j ] ; 
	}
	for( i = *div ; i < *n ; i++ ) 
	{
		for( j = 0 ; j < *m ; j++ ) 
			out2[ i - *div + (*n - *div) * j ] = in[ i + *n * j ] ; 
	}
}

// Translates a matrix of ints to doubles 
// in : int matrix 
// n : rows 
// m : cols 
// out : a double matrix with rows X cols entries 
void intToDoubleMat ( int *in , int *n , int *m , double *out ) 
{
	int i , j ; 
	for( i = 0 ; i < *n ; i++ )
	{
		for( j = 0 ; j < *m ; j++ )
			out[ i + *n * j ] = (double) in[ i + *n * j ] ; 
	}
}

// Calculates the minor of a matrix 
// in : a column-indexed n X n matrix 
// n : pointer to an int storing the matrix dimension 
// i : the row to be removed 
// j : the column to be removed 
// out : an address for a single double to be stored 
void minorVal ( double *in , int *n , int *i , int *j , double *out ) 
{
	if( *n < 3 )
	{
		printf( "ERROR: minor of small matrix\n" ) ; 
		return ; 
	}
	
	double *tmp = (double*) malloc( (*n-1)*(*n-1)*sizeof(double) ) ; 
	int ii , jj ; 
	int m = *n - 1 ;  
	for( ii = 0 ; ii < *n ; ii++ ) 
	{
		for( jj = 0 ; jj < *n ; jj++ ) 
		{
			if( ii < *i && jj < *j ) 
				tmp[ ii + m * jj ] = in[ ii + *n * jj ] ; 
			if( ii > *i && jj < *j ) 
				tmp[ ii-1 + m * jj ] = in[ ii + *n * jj ] ; 
			if( ii < *i && jj > *j ) 
				tmp[ ii + m * ( jj - 1 ) ] = in[ ii + *n * jj ] ; 
			if( ii > *i && jj > *j ) 
				tmp[ ii - 1 + m * ( jj - 1) ] = in[ ii + *n * jj ] ; 
		}
	}
	det( tmp , &m , out ) ; 
	free( tmp ) ; 
} 

// Calculates the determinant of matrix 
// in : a column-indexed n X n matrix 
// out : an address for a single double 
void det ( double *in , int *n , double *out ) 
{
	if( *n == 0 ) 
	{
		out = NULL ; 
		return ; 
	}
	if( *n == 1 ) 
	{
		*out = *in ; 
		return ; 
	}
	if( *n == 2 ) 
	{
		*out = in[0] * in[3] - in[1] * in[2] ; 
		return ; 
	}
	int i ; 
	*out = 0.0 ; 
	double sign = 1.0 ; 
	double min ; 
	int zero = 0 ; 
	for( i = 0 ; i < *n ; i++ ) 
	{
		minorVal( in , n , &i , &zero , &min ) ; 
		*out += sign * in[ i ] * min ; 
		sign = -sign ; 
	}
} 

// Calculates the adjoint matrix 
// in : a column-indexed n X n matrix 
// out : enough space for an n X n matrix 
void adj ( double *in , int *n , double *out ) 
{
	int i , j ; 
	double min ; 
	for( i = 0 ; i < *n ; i++ ) 
	{
		for( j = 0 ; j < *n ; j++ ) 
		{
			minorVal( in , n , &i , &j , &min ) ; 
			out[ j + *n * i ] = pow( -1.0 , (double) (i + j + 2) ) * min ; 
		}
	}
} 

// NOT NUMERICALLY STABLE 
// Calculates the inverse of a matrix 
// in : a column-indexed n X n matrix 
// out : enough space for an n X n matrix 
void inv ( double *in , int *n , double *out ) 
{
	if( *n == 0 ) 
	{
		out = NULL ; 
		return ; 
	}
	if( *n == 1 ) 
	{
		*out = *in ; 
		return ; 
	}
	double c ; 
	if( *n == 2 ) 
	{
		c = in[0] * in[3] - in[1] * in[2] ; 
		out[0] = in[3] / c ; 
		out[1] = -in[1] / c ; 
		out[2] = -in[2] / c ; 
		out[3] = in[0] / c ; 
		return ; 
	}
	det( in , n , &c ) ; 
	adj( in , n , out ) ; 
	int i , m ; 
	m = (*n) * (*n) ; 
	for( i = 0 ; i < m ; i++ ) 
	{
		out[i] = out[i] / c ; 
	}
}

struct matProdArg 
{
	double *a ; 
	double *b ; 
	int *m ; 
	int *n ; 
	int *p ; 
	double *out ; 
	int start ; // job id 
	int stop ; 
}; 

void* matProdPosix ( void * targ ) 
{
	struct matProdArg *arg = (struct matProdArg*) targ ; 
	int i , j , k , l ; 
	if( *(arg->m) > 1 && *(arg->p) > 1 ) 
	{
		for( l = arg->start ; l < arg->stop ; l++ ) // job id 
		{
			i = floor( ((double) l) / ((double) *(arg->m)) ) ; 
			j = l % *(arg->m) ; 
			arg->out[ i + *(arg->m) * j ] = 0.0 ; 
			for( k = 0 ; k < *(arg->m) ; k++ ) 
				arg->out[ i + *(arg->m) * j ] += arg->a[ i + *(arg->m) * k ] * arg->b[ k + *(arg->n) * j ] ; 
		}
	}
	
	if( *(arg->m) > 1 && *(arg->p) == 1 )
	{
		j = 0 ; 
		for( i = arg->start ; i < arg->stop ; i++ ) 
		{
			arg->out[ i ] = 0.0 ; 
			for( k = 0 ; k < *(arg->m) ; k++ ) 
				arg->out[ i ] += arg->a[ i + *(arg->m) * k ] * arg->b[ k ] ; 
		}
	}
	
	if( *(arg->m) == 1 && *(arg->p) > 1 )
	{
		i = 0 ; 
		for( j = arg->start ; j < arg->stop ; j++ ) 
		{
			arg->out[ j ] = 0.0 ; 
			for( k = 0 ; k < *(arg->p) ; k++ ) 
				arg->out[ j ] += arg->a[ k ] * arg->b[ k + *(arg->n) * j ] ; 
		}
	}
	
	if( *(arg->m) == 1 && *(arg->p) == 1 ) 
	{
		arg->out[0] = 0.0 ; 
		for( k = 0 ; k < *(arg->m) ; k++ ) 
			arg->out[0] += arg->a[k] * arg->b[ k ] ; 
	}
	
	return NULL ; 
}

// Calculates the product of matrices a & b in parallel  
// a : a column-indexed m X n matrix 
// b : a column-indexed n X p matrix 
// out : enough space for an m X p matrix 
// threads : the number of threads to utilize 
void matProdP ( double *a , double *b , int *m , int *n , int *p , double *out , int *threads ) 
{
	// allocate work 
	int totalWork = *m * *p ; 
	int thrN = *threads ; 
	if( thrN > totalWork ) 
		thrN = totalWork ; 
	pthread_t *thr = (pthread_t*) malloc( thrN * sizeof(pthread_t) ) ; 
	struct matProdArg *parg = (struct matProdArg*) malloc( thrN * sizeof(struct matProdArg) ) ; 
	int i ; 
	for( i = 0 ; i < thrN ; i++ ) 
	{
		parg[i].a = a ; 
		parg[i].b = b ; 
		parg[i].m = m ; 
		parg[i].n = n ; 
		parg[i].p = p ; 
		parg[i].out = out ; 
		parg[i].start = floor( ((double) i * totalWork) / ((double) thrN) ) ; 
		parg[i].stop = floor( ((double) (i+1) * totalWork) / ((double) thrN) ) ; 
		if( i + 1 == thrN ) 
			parg[i].stop = totalWork ; 
	}
	
	// run threads 
	int err ; 
	for( i = 0 ; i < thrN ; i++ ) 
	{
		err = pthread_create( &thr[i] , NULL , &matProdPosix , (void*) &parg[i] ) ; 
		if( err ) 
		{
			printf( "ERROR, POSIX: %i\n" , err ) ; 
			return ; 
		}
	}
	
	// join threads 
	void *status ; 
	for( i = 0 ; i < thrN ; i++ ) 
	{
		err = pthread_join( thr[i] , &status ) ; 
		if( err ) 
		{
			printf( "ERROR, POSIX: %i\n" , err ) ;  
			return ; 
		}
	}
	
	// job's done 
}

// Calculates the product of matrices a & b 
// a : a column-indexed m X n matrix 
// b : a column-indexed n X p matrix 
// out : enough space for an m X p matrix 
void matProd ( double *a , double *b , int *m , int *n , int *p , double *out )
{
        int i , j , k ;
        for( i = 0 ; i < *m ; i++ )
        {
                for( j = 0 ; j < *p ; j++ )
                {
                        out[ i + *m * j ] = 0.0 ;
                        for( k = 0 ; k < *n ; k++ )
                                out[ i + *m * j ] += a[ i + *m * k ] * b[ k + *n * j ] ;
                }
        }
}

// Calculates the product of aT b
// a : a column-indexed n X m matrix 
// b : a column-indexed n X p matrix 
// out : enough space for an m X p matrix 
void matTrProd ( double *a , double *b , int *m , int *n , int *p , double *out ) 
{
	int i , j , k ; 
	for( i = 0 ; i < *m ; i++ ) 
	{
		for( j = 0 ; j < *p ; j++ ) 
		{
			out[ i + *m * j ] = 0.0 ; 
			for( k = 0 ; k < *m ; k++ ) 
				out[ i + *m * j ] += a[ k + *n * i ] * b[ k + *n * j ] ; 
		}
	}
}

// Calculates eigenvalues via the QR algorithm 
// in : a column-indexed n X n matrix 
// eps : convergence parameter, choose somthing small 
// maxIter : maximum number of iterations allowed 
// val : space for n doubles, the eigen values 
// vec : space for an n X n matrix, the eigen vectors 
void qrEig ( double *in , int *n , double *eps , int *minIter , int *maxIter , double *val , double *vec ) 
{
	double err = *eps + 1.0 ; 
	double *prevVal = (double*) malloc( *n * sizeof(double) ) ; 
	double *prevVec = (double*) malloc( *n * *n * sizeof(double) ) ; 
	double *aMat = (double*) malloc( *n * *n * sizeof(double) ) ; 
	double *qMat = (double*) malloc( *n * *n * sizeof(double) ) ; 
	double *rMat = (double*) malloc( *n * *n * sizeof(double) ) ; 
	double *uMat = (double*) malloc( *n * *n * sizeof(double) ) ; 
	int k = -1 ; 
	int i , j ; 
	for( i = 0 ; i < *n ; i++ ) 
	{
		for( j = 0 ; j < *n ; j++ ) 
		{
			aMat[ i + *n * j ] = in[ i + *n * j ] ; 
			if( i == j )
				uMat[ i + *n * i ] = 1.0 ; 
			else
				uMat[ i + *n * j ] = 0.0 ; 
		}
	}
	while( (err > *eps || k < *minIter ) && k < *maxIter ) 
	{
		// Copy previous eigen values 
		for( i = 0 ; i < *n ; i++ ) 
			prevVal[i] = val[i] ; 
		
		// QR Decomp 
		qrDecomp( aMat , n , qMat , rMat ) ; 
		
		// update Schur decomposition   
		matProd( rMat , qMat , n , n , n , aMat ) ; 
		
		// update vectors 
		matProd( uMat , qMat , n , n , n , vec ) ; 
		for( i = 0 ; i < *n ; i++ ) 
		{
			for( j = 0 ; j < *n ; j++ ) 
			{
				prevVec[ i + *n * j ] = uMat[ i + *n * j ] ; 
				uMat[ i + *n * j ] = vec[ i + *n * j ] ; 
			}
		}
		
		// update values 
		for( i = 0 ; i < *n ; i++ )  
			val[i] = aMat[ i + *n * i ] ; 
		
		// update iterations 
		k = k + 1 ; 
		if( k > 0 ) 
		{
			err = 0.0 ; 
			for( i = 0 ; i < *n ; i++ ) 
			{
				err += fabs( prevVal[i] - val[i] ) ; 
				for( j = 0 ; j < *n ; j++ )
					err += fabs( prevVec[ i + *n * j ] - uMat[ i + *n * j ] ) ; 
			}
		}
		
	}
	*maxIter = k ; 
	for( i = 0 ; i < *n ; i++ )
	{
		for( j = 0 ; j < *n ; j++ ) 
			vec[ i + *n * j ] = uMat[ i + *n * j ] ; 
	}
	
	free( prevVal ) ; 
	free( aMat ) ; 
	free( qMat ) ; 
	free( uMat ) ; 
} 

// Calculates the QR decomposition 
// in : a column-indexed n X n matrix 
// q : space for an n X n matrix 
// r : space for an n X n matrix 
void qrDecomp ( double *in , int *n , double *q , double *r ) 
{
	gramSchmidt( in , n , q ) ; 
	int i , j , k ; 
	for( i = 0 ; i < *n ; i++ ) 
	{
		for( j = 0 ; j < *n ; j++ ) 
		{
			r[ i + *n * j ] = 0.0 ; 
			if( j >= i ) 
			{
				for( k = 0 ; k < *n ; k++ ) 
				{
					r[ i + *n * j ] += q[ k + *n * i ] * in[ k + *n * j ] ; 
				}
			}
		}
	}
} 

// Calculates an orthogonal basis via the Gram-Schmidt algorithm 
// in : a column-indexed n X n matrix 
// out : space for an n X n matrix 
void gramSchmidt ( double *in , int *n , double *out ) 
{
	int i , j , k ; 
	double *tmp = (double*) malloc( *n * sizeof(double) ) ; 
	double total ; 
	for( i = 0 ; i < *n ; i++ ) // for each basis vector 
	{
		// copy initial vector 
		for( j = 0 ; j < *n ; j++ ) 
			out[ j + *n * i ] = in[ j + *n * i ] ; 
		
		for( j = 0 ; j < i ; j++ ) // for each basis vector before this one 
		{
			// perject this vector onto the basis vector 
			proj( &in[ *n * i ] , &out[ *n * j ] , n , tmp ) ; 
			
			// subtract the projection from this vector 
			for( k = 0 ; k < *n ; k++ ) 
				out[ k + *n * i ] -= tmp[k] ;  
		}
		
		// standardize this vector to make it a unit-length basis vector 
		total = 0.0 ; 
		for( j = 0 ; j < *n ; j++ ) 
			total += out[ j + *n * i ] * out[ j + *n * i ] ; 
		total = sqrt(total) ; 
		for( j = 0 ; j < *n ; j++ ) 
			out[j + *n * i ] = out[ j + *n * i ] / total ; 
	} 
	free(tmp) ; 
}

// Calculates the projection of vector x onto y 
void proj ( double *x , double *y , int *n , double *out ) 
{
	double tmp1 = 0.0 ; 
	double tmp2 = 0.0 ; 
	int i ; 
	for( i = 0 ; i < *n ; i++ ) 
	{
		tmp1 += x[i] * y[i] ; 
		tmp2 += y[i] * y[i] ; 
	}
	for( i = 0 ; i < *n ; i++ ) 
		out[i] = tmp1 * y[i] / tmp2 ; 
}

// Appends two matrices 
// top : the matrix desired on the top 
// bot : the matrix desired on the bottom 
// rows1 : the number of rows in top 
// rows2 : the number of rows in bottom 
// cols : the number of columns in both matrices 
// out : enough space for (rows1 + rows2) * cols doubles 
void appendRows ( double *top , double *bot , int *rows1 , int *rows2 , int *cols , double *out ) 
{
	int i , j ; 
	for( i = 0 ; i < *rows1 ; i++ ) 
	{
		for( j = 0 ; j < *cols ; j++ ) 
			out[ i + (*rows1 + *rows2) * j ] = top[ i + *rows1 * j ] ; 
	}
	for( i = 0 ; i < *rows2 ; i++ ) 
	{
		for( j = 0 ; j < *cols ; j++ ) 
			out[ i + *rows1 + (*rows1 + *rows2) * j ] = bot[ i + *rows2 * j ] ; 
	}
}

// Calculates the Cholesky decomposition
// in : a column-major positive-definite n X n matrix 
// n : > 1 
// out : enough space for an n X n matrix 
void chol ( double *in , int *n , double *out ) 
{
	int i , j , k ; 
	double tmp ; 
	for( i = 0 ; i < *n ; i++ ) 
	{
		for( j = 0 ; j < *n ; j++ ) 
		{
			if( i < j )
			{
				out[ i + *n * j ] = 0.0 ; 
			}
			else
			{
				tmp = in[i + *n * j] ; 
				if( i == j ) 
				{
					for( k = 0 ; k < j ; k++ )
					{
						tmp -= out[ i + *n * k ] * out[ i + *n * k ] ; 
					}
					out[ i + *n * j ] = sqrt(tmp) ; 
				}
				else
				{
					for( k = 0 ; k < j ; k++ )
					{
						tmp -= out[i + *n * k] * out[j + *n * k] ; 
					} 
					out[i + *n * j] = tmp / out[j + *n * j] ; 
				}
			}
		}
	}
}

void detPsd ( double *in , int *n , double *out ) 
{
	double *tmp = (double*) malloc( *n * *n * sizeof(double) ) ; 
	chol( in , n , tmp ) ; 
	*out = 0.0 ; 
	int i ; 
	for( i = 0 ; i < *n ; i++ )
		*out += tmp[ i + *n * i ] ;
	*out = *out * *out ;  
	free(tmp) ; 
}

// calculates the inverse matrix of a positive-symmetric-definate matrix 
// in : a column-major order n X n PSD matrix to be inverted 
// out : space for n X n doubles 
// eps : convergence term 
// threads : number of available POSIX threads 
void invPsd( double *in , int *n , double *out , double *eps , int *threads )  
{
	double *val = (double*) malloc( *n * sizeof(double) ) ; 
	double *vec = (double*) malloc( *n * *n * sizeof(double) ) ; 
	double *diag = (double*) malloc( *n * *n * sizeof(double) ) ; 
	double *tmp = (double*) malloc( *n * *n * sizeof(double) ) ; 
	// powerIteration( in , n , n , eps , val , vec , threads ) ; 
	double eps2 = *eps * 0.0001 ; 
	psdEig ( in , n , &eps2 , vec , val ) ; 
	
	int i , j ; 
	for( i = 0 ; i < *n ; i++ ) 
	{
		for( j = 0 ; j < *n ; j++ ) 
		{
			if( i == j ) 
				diag[ i + *n * j ] = 1.0 / val[i] ; 
			else
				diag[ i + *n * j ] = 0.0 ; 
		}
	}
	matProd ( vec , diag , n , n , n , tmp ) ; // TODO debug matProdP and upgrade to parallel  
	transpose ( vec , n , n , diag ) ; // use diag as working space 
	matProd ( tmp , diag , n , n , n , out ) ; 
	
	free( val ) ; 
	free( vec ) ; 
	free( diag ) ; 
	free( tmp ) ; 
}

// converts a symmetric matrix to a tri-diagonal matrix through a series of similar transforms 
// mat : a column-major n X n matrix 
// tri : enough space for n*n doubles 
// q : enough space for n*n doubles 
// void symmToTri ( double *mat , int *n , double *tri , double *q , int *threads ) 
void symmToTri ( double *mat , int *n , double *tri ) 
{
	double *u = (double*) malloc( *n * sizeof(double) ) ; 
	double *v = (double*) malloc( *n * sizeof(double) ) ; 
//	double *p = (double*) malloc( *n * *n * sizeof(double) ) ; 
//	double *tmpMat = (double*) malloc( *n * *n * sizeof(double) ) ; 
//	double *tmpPtr ; 
//	double *firstQ = q ; // I utilize a pointer-swap for speed-up, but must return on the correct address 
//	double *firstTmpMat = tmpMat ; 
	int i , j , k ; 
	for( i = 0 ; i < *n ; i++ ) 
	{
		for( j = 0 ; j < *n ; j++ ) 
		{
			tri[ i + *n * j ] = mat[ i + *n * j ] ; 
			/* 
			if( i == j ) 
				q[ i + *n * j ] = 1.0 ; 
			else
				q[ i + *n * j ] = 0.0 ; 
			*/
		}
	}
	
	double tmp ; 
	for( k = 1 ; k < *n - 1 ; k++ ) 
	{
		// design u, of length n - k
		tmp = 0.0 ;  
		for( i = 0 ; i < *n - k ; i++ ) 
		{
			u[ i ] = tri[ k + i + *n * (k-1) ] ; 
			tmp += u[i] * u[i] ; 
		}
		u[0] -= copysign(1.0,u[0]) * sqrt(tmp) ; 
		tmp = 0.0 ; 
		for( i = 0 ; i < *n - k ; i++ ) 
			tmp += u[i] * u[i] ; 
		for( i = 0 ; i < *n - k ; i++ ) 
			u[i] = u[i] / sqrt(tmp) ; 
		
		// design p for later calculation of q 
		/*
		for( i = 0 ; i < *n ; i++ ) 
		{
			for( j = 0 ; j < *n ; j++ ) 
			{
				if( i == j ) 
					p[ i + *n * j ] = 1.0 ; 
				else
					p[ i + *n * j ] = 0.0 ; 
				if( i  >= k && j >= k ) 
					p[ i + *n * j ] -= 2.0 * u[ i - k ] * u[ j - k ] ; 
			}
		}
		*/
		
		// design v 
		tmp = 0.0 ; 
		for( i = 0 ; i < *n - k ; i++ ) 
		{
			for( j = 0 ; j < *n - k ; j++ ) 
				tmp += tri[ i + k + *n * (j + k) ] * u[i] * u[j] ; // quadratic form 
		}
		for( i = 0 ; i < *n ; i++ ) 
			v[i] = 0.0 ; 
		for( i = 0 ; i < *n ; i++ ) 
		{
			for( j = k ; j < *n ; j++ ) 
				v[i] += tri[ i + *n * j ] * u[ j - k ] ; 
			if( i < *n - k )
				v[i+k] -= u[i] * tmp ; 
			v[i] = 2.0 * v[i] ; 
		}
		
		// advance tri closer to a tri-diagonal 
		for( i = 0 ; i < *n ; i++ ) 
		{
			for( j = 0 ; j < *n ; j++ ) 
			{
				if( (j == k-1 && i > k ) || (i == k-1 && j > k ) ) 
				{
					tri[ i + *n * j ] = 0.0 ; 
				}
				else
				{
					if( i >= k && j < k ) 
						tri[ i + *n * j ] -= u[i-k] * v[j] ; 
					if( i < k && j >= k ) 
						tri[ i + *n * j ] -= v[i] * u[j-k] ; 
					if( i >= k && j >= k ) 
						tri[ i + *n * j ] -= u[i-k] * v[j] + v[i] * u[j-k] ; 
				}
			}
		}
		
		/*
		// advance q closer to its final form 
		matProdP ( q , p , n , n , n , tmpMat , threads ) ; 
		tmpPtr = q ; 
		q = tmpMat ; // pointer swap 
		tmpMat = tmpPtr ; 
		*/
		
	}
	
	// ensure that q points to the correct address after pointer swaps 
//	if( firstQ != q ) 
//		memcpy( (void*) firstQ , (void*) q , *n * *n * sizeof(double) ) ; 
	
	free( u ) ; 
	free( v ) ; 
//	free( p ) ; 
//	free( firstTmpMat ) ; 
}

// iteration function for power iteration method 
// tri : an n X n tri-diagonal matrix in column major order 
// prev : a unit-length n X 1 vector 
// post : space for a unit-length n X 1 vector 
// lam : estimate for the eigen value 
void powIterTridiag( double *tri , int *n , double *prev , double *post , double *lam ) 
{
	// first & last entries, unscaled  
	post[0] = prev[0] * tri[0] + prev[1] * tri[ *n ] ; 
	*lam = post[0] * post[0] ; 
	post[*n - 1] = prev[*n - 2] * tri[ *n-1 + *n *(*n - 2) ] + prev[*n - 1] * tri[ *n - 1 + *n * (*n - 1) ] ;  
	*lam += post[*n - 1] * post[*n - 1] ; 
	// mid entries, unscaled 
	int i ; 
	for( i = 1 ; i < *n - 1 ; i++ ) 
	{
		post[i] = prev[i-1] * tri[ i + *n * (i-1) ] + prev[i] * tri[ i + *n * i ] + prev[i+1] * tri[ i + *n * (i+1) ] ; 
		*lam += post[i] * post[i] ; 
	}
	
	// rescale entries  
	*lam = sqrt( *lam ) ; 
	for( i = 0 ; i < *n ; i++ ) 
		post[i] = post[i] / *lam ; 
}

void getEigenPowerItr( double *tri , int *n , double *eps , double *val , double *vec , int *threads ) 
{
	// random initial vector 
	int i ; 
	*val = 0.0 ; 
	for( i = 0 ; i < *n ; i++ ) 
	{
		vec[i] =  ((double) ( rand() % 100000 ) / 100000.0 ) ; // TODO check zero case  
		*val += vec[i] * vec[i] ; 
	}
	// scale to unit length 
	*val = sqrt( *val ) ; 
	for( i = 0 ; i < *n ; i++ ) 
		vec[i] = vec[i] / *val ; 
// printf( "DEBUG, initial vec: " ) ; for( i = 0 ; i < *n ; i++ ){ printf( "%f " , vec[i] ) ; } ; printf( "\n" ) ; 
	
	// begin power iteration 
	double *tmpVec = (double*) malloc( *n * sizeof(double) ) ; 
	double err = *eps + 1.0 ; 
	double prev ; 
	int one = 1.0 ; 
	while( err > *eps ) 
	{
		prev = *val ; 
		
//		*val = 0.0 ; 
//		powIterTridiag( tri , n , vec , tmpVec , val ) ; // pointer trick 
//		powIterTridiag( tri , n , tmpVec , vec , val ) ; 
		
		matProd( tri , vec , n , n , &one , tmpVec ) ; 
		*val = 0.0 ; 
		for( i = 0 ; i < *n ; i++ ) 
			*val += tmpVec[i] * tmpVec[i] ; 
		*val = sqrt( *val ) ; 
		for( i = 0 ; i < *n ; i++ ) 
			tmpVec[i] = tmpVec[i] / *val ; 
// printf( "DEBUG, tmpVec: " ) ; for( i = 0 ; i < *n ; i++ ){ printf( "%f " , tmpVec[i] ) ; } ; printf( "\n" ) ;
		
		matProd( tri , tmpVec , n , n , &one , vec ) ; // pointer trick 
		*val = 0.0 ; 
		for( i = 0 ; i < *n ; i++ ) 
			*val += vec[i] * vec[i] ; 
		*val = sqrt( *val ) ; 
		for( i = 0 ; i < *n ; i++ ) 
			vec[i] = vec[i] / *val ; 
		
		err = fabs( *val - prev ) ; 
// printf( "DEBUG, val: %f, err: %f\n" , *val , err ) ; 
// printf( "DEBUG, vec: " ) ; for( i = 0 ; i < *n ; i++ ){ printf( "%f " , vec[i] ) ; } ; printf( "\n" ) ; 
	}
	
	free( tmpVec ) ; 
}

// removes all co-linearity of colspace(x) with y 
// x : n X n matrix 
// y : an n X 1 vector 
// out : space for an n X n matrix 
void removeDimension( double *x , int *n , double *y , double *out )
{
        int i , j ;
        double *tmp = (double*) malloc( *n * sizeof(double) ) ;
        for( i = 0 ; i < *n ; i++ )
        {
                proj( &x[ *n * i ] , y , n , tmp ) ;
                for( j = 0 ; j < *n ; j++ )
                        out[ j + *n * i ] = x[ j + *n * i] - tmp[j] ;
        }
        free( tmp ) ;
}

// calculates eigen vecs & vals via power iteration on a symmetric matrix 
// accurracy may suffer slightly in exchange for more speed 
// mat : a symmetric n X n matrix 
// m : the number of desired eigen vectors 
// eps : tolerance, choose something small like e-14 
// val : output space, enough space for m eigen values 
// vec : output space, enough space for n X m doubles 
// threads : number of POSIX threads desired for the calculation 
void powerIteration( double *mat , int *n , int *m , double *eps , double *val , double *vec , int *threads ) 
{
	double *tri = (double*) malloc( *n * *n * sizeof(double) ) ; 
	double *tmpMat = (double*) malloc( *n * *n * sizeof(double) ) ; 
	double *adjMat = (double*) malloc( *n * *n * sizeof(double) ) ; 
	double *tmpVec = (double*) malloc( *n * sizeof(double) ) ; 
	
	memcpy( tmpMat , mat , *n * *n * sizeof(double) ) ; 
	
	int i , j , k ; 
	for( i = 0 ; i < *m ; i++ ) 
	{
		/* 
		// convert the previous matrix to a similar tri-diagonal matrix 
		if( i == 0 ) 
			symmToTri( mat , n , tri ) ; 
		else
			symmToTri( tmpMat , n , tri ) ; 
		*/ 
		
		// find the largest eigenvalue and its vector 
		getEigenPowerItr( tmpMat , n , eps , &val[i] , &vec[ *n * i ] , threads ) ;  
		
		// remove the eigenvector from the span if another is wanted 
		if( i < *m - 1 ) 
		{
			removeDimension( tmpMat , n , &vec[*n * i] , adjMat ) ;  
			memcpy( tmpMat , adjMat , *n * *n * sizeof(double) ) ;  
		}
	}
	
	/* 
	// for each eigenvalue, calculate its vector 
	// this should be done in parallel TODO 
	for( i = 0 ; i < *m ; i++ ) 
	{
		// adjust the matrix for the eigenvalue 
		memcpy( adjMat , mat , *n * *n * sizeof(double) ) ; 
		for( j = 0 ; j < *n ; j++ ) 
			adjMat[ j + *n * j ] -= val[i] ; 
		
		// calculate the cholesky decomp 
		chol( adjMat , n , tmpMat ) ; 
		
		// WLOG, 
	}
	*/
	
	free( tri ) ; 
	free( tmpMat ) ; 
	free( adjMat ) ; 
	free( tmpVec ) ; 
}

double min( double x , double y ) 
{
	if( x < y ) 
		return x ; 
	return y ; 
}

double max( double x , double y ) 
{
	if( x < y ) 
		return y ; 
	return x ; 
}

double hypot( double x , double y ) 
{
	x = fabs(x) ; 
	y = fabs(y) ; 
	double t = min(x,y) ; 
	x = max(x,y) ; 
	t = t / x ; 
	return x * sqrt( 1.0 + t * t ) ; 
}

void givens( double *x , double *y , double *c , double *s ) 
{
	double r = hypot( *x , *y ) ; 
	*c = *x / r ; 
	*s = - *y / r ; 
}

// provides rotation vectors for a Householder rotation 
// x : a col-major order PSD matrix with rows rows 
// col : the column to rotate 
// u , v : output vectors of lenth rows 
void householder ( double *x , int *rows , int *col , double *u , double *v ) 
{
	int i , j ; 
	double ttl = 0.0 ; 
	for( i = 0 ; i < *col + 1 ; i++ ) 
		u[i] = 0.0 ; 
	for( i = *col + 1 ; i < *rows ; i++ ) 
	{
		u[i] = x[ i + *rows * *col ] ; 
		ttl += x[ i + *rows * *col ] * x[ i + *rows * *col ] ; 
	}
	ttl = sqrt(ttl) ; 
	u[ *col + 1 ] += copysign( ttl , u[ *col + 1 ] ) ; 
	ttl = 0.0 ; 
	for( i = *col + 1 ; i < *rows ; i++ ) 
		ttl += u[i] * u[i] ; 
	ttl = sqrt(ttl) ; 
	for( i = 0 ; i < *rows ; i++ ) 
		u[i] = u[i] / ttl ; // u complete  
	if( v == NULL ) 
		return ; 
	// calculate v , for PSD x only , v = 2 * x %*% u - 2 * u * ( t(u) %*% x %*% u )  
	ttl = 0.0 ; // this will hold a quadratic form 
	for( i = *col + 1 ; i < *rows ; i++ ) // calculate avoiding multiplying zeros 
	{
		for( j = *col + 1 ; j < *rows ; j++ ) 
			ttl += x[ i + *rows * j ] * u[i] * u[j] ; 
	}
	for( i = 0 ; i < *rows ; i++ ) 
		v[i] = -2.0 * u[i] * ttl ; 
	for( i = 0 ; i < *rows ; i ++ ) 
	{
		ttl = 0.0 ; 
		for( j = *col + 1 ; j < *rows ; j++ ) 
			ttl += x[ i + *rows * j ] * u[j] ; 
		ttl *= 2.0 ; 
		v[i] += ttl ; 
	} 
}
// Eigen calculator for PSD matrices 
// x : column-major order n X n PSD matrix 
// eps : convergence term 
// q : space for n X n doubles for eigenvectors, if null the algorithm will be sped up but no eigenvectors will be returned 
// vals : space for n doubles for eigenvalues 
void psdEig ( double *x , int *n , double *eps , double *q , double *vals ) 
{
	double *y = (double*) malloc( *n * *n * sizeof(double) ) ; 
	memcpy( y , x , *n * *n * sizeof(double) ) ; 
	int *idx = (int*) malloc( (*n) * sizeof(int) ) ; 
          
	double *qq = NULL ; 
	int i , j , k , l ; 
	if( q != NULL ) 
	{
		qq = (double*) malloc( *n * *n * sizeof(double) ) ; 
		for( i = 0 ; i < *n ; i++ ) 
		{
			for( j = 0 ; j < *n ; j++ ) 
			{
				if( i == j ) 
					q[ i + *n * j ] = 1.0 ; 
				else
					q[ i + *n * j ] = 0.0 ; 
			}
		}
	}
	
	// householder rotations 
	double *u = (double*) malloc( *n * sizeof(double) ) ; 
	double *v = (double*) malloc( *n * sizeof(double) ) ; 
	for( i = 0 ; i < *n - 1 ; i++ ) 
	{
		householder ( y , n , &i , u , v ) ; 
		for( j = i+1 ; j < *n ; j++ ) 
		{
			for( k = 0 ; k < *n ; k++ ) 
				y[ j + *n * k ] -= u[j] * v[k] ; 
		} 
		for( j = 0 ; j < *n ; j++ ) 
		{
			for( k = i+1 ; k < *n ; k++ ) 
				y[ j + *n * k ] -= v[j] * u[k] ; 
		}
		if( q != NULL ) 
		{
			for( j = 0 ; j < *n ; j++ ) 
			{
				for( k = i+1 ; k < *n ; k++ ) 
				{ 
					qq[ j + *n * k ] = 0.0 ; 
					for( l = i+1 ; l < *n ; l++ ) 
						qq[ j + *n * k ] += q[ j + *n * l ] * u[l] * u[k] ; 
					qq[ j + *n * k ] *= 2.0 ; 
				}
			}
			for( j = 0 ; j < *n ; j++ ) 
			{
				for( k = i+1 ; k < *n ; k++ ) 
					q[ j + *n * k ] -= qq[ j + *n * k ] ; 
			}
		}
	}
	
	// symmetric tri-diagonal qr algorithm with implicit Wilkinson shifts 
	double err , d , s , c , w , z , xx , yy , t ; 
	for( i = *n ; i > 1 ; i-- ) 
	{
		err = 1.0 + *eps ; 
		while( err > *eps ) 
		{
// fprintf( stderr , "DEBUG, err: %e, eps: %e, i: %i\n" , err , *eps , i ) ; 
			d = ( y[i-2 + *n * (i-2)] - y[i-1 + *n * (i-1)] ) / 2.0 ; 
			if( d == 0.0 ) 
				s = y[i-1 + *n * (i-1)] - fabs( y[ i-2 + *n * (i-1) ] ) ; 
			else
				s = y[i-1 + *n * (i-1)] - y[ i-2 + *n * (i-1) ] * y[ i-2 + *n * (i-1) ] / ( d + copysign( sqrt( d*d + y[ i-2 + *n * (i-1) ] * y[ i-2 + *n * (i-1) ] ) , d ) ) ; 
			xx = y[0] - s ; 
			yy = y[ *n ] ; 
			for( j = 0 ; j < i-1 ; j++ ) 
			{
				if( i > 1 ) 
					givens ( &xx , &yy , &c , &s ) ; 
				else 
				{
					t = 0.5 * atan( 2.0 * y[*n] / ( y[ 1 + *n ] - y[0] ) ) ; 
					c = cos(t) ; 
					s = sin(t) ; 
				}
				w = c * xx - s * yy ; 
				d = y[ j + *n * j ] - y[ j+1 + *n * (j+1) ] ; 
				z = s * ( 2.0 * c * y[ j + *n * (j+1) ] + d * s ) ; 
				y[ j + *n * j ] = y[ j + *n * j ] - z ; 
				y[ j+1 + *n * (j+1) ] = y[ j+1 + *n * (j+1) ] + z ; 
				y[ j + *n * (j+1) ] = d * c * s + (c*c - s*s) * y[ j + *n * (j+1) ] ; 
				y[ j+1 + *n * j ] = y[ j + *n * (j+1) ] ; 
				xx = y[ j + *n * (j+1) ] ; 
				if( j > 0 ) 
				{
					y[ j-1 + *n * j ] = w ; 
					y[ j + *n * (j-1) ] = w ; 
				}
				if( j < i-1 ) 
				{
					yy = -s * y[ j+1 + *n * (j+2) ] ; 
					y[ j+1 + *n * (j+2) ] = c * y[ j+1 + *n * (j+2) ] ; 
					y[ j+2 + *n * (j+1) ] = y[ j+1 + *n * (j+2) ] ; 
				} 
				
				// update vectors 
				if( q != NULL ) 
				{
					for( k = 0 ; k < *n ; k++ ) 
					{
						qq[ k ] = c * q[ k + *n * j ] - s * q[ k + *n * (j+1) ] ; 
						qq[ k + *n ] = s * q[ k + *n * j ] + c * q[ k + *n * (j+1) ] ;  
					}
					for( k = 0 ; k < *n ; k++ ) 
					{
						q[ k + *n * j ] = qq[ k ] ; 
						q[ k + *n * (j+1) ] = qq[ k + *n ] ; 
					}
				}
				
				// calculate convergence errors 
				// err = abs( y[ i-2 + *n * (i-1) ] ) / ( abs( y[ i-2 + *n * (i-2) ] ) + abs( y[ i-1 + *n * (i-1) ] ) ) ; 
			}
			// calculate convergence errors 
			err = fabs( y[ i-2 + *n * (i-1) ] ) / ( fabs( y[ i-2 + *n * (i-2) ] ) + fabs( y[ i-1 + *n * (i-1) ] ) ) ; 
		}
	}
	
	// eigenvalues are stored on diagonal 
	for( i = 0 ; i < *n ; i++ ) 
	{
		idx[i] = i ; 
		vals[i] = y[ i + *n * i ] ; 
		y[i] = y[ i + *n * i ] ; 
		if( q != NULL ) 
		{
			for( j = 0 ; j < *n ; j++ ) 
				qq[ i + *n * j ] = q[ i + *n * j ] ; 
		}
	}
	
	// sort eigenvalues, TODO implement as pivoting  
	quickSort( vals , n , idx ) ; 
	for( i = 0 ; i < *n ; i++ ) 
	{
		vals[*n - i - 1] = y[ idx[i] ] ; 
		if( q != NULL ) 
		{
			for( j = 0 ; j < *n ; j++ ) 
				q[ j + *n * (*n - i - 1) ] = qq[ j + *n * idx[i] ] ; 
		}
	}
	
//fprintf( stderr , "DEBUG 1 u is from %p to %p\n" , u , u + *n * sizeof(double) ) ; 
	free( u ) ; 
//fprintf( stderr , "DEBUG 2 v is from %p to %p\n" , v , v + *n * sizeof(double) ) ;
	free( v ) ; 
//fprintf( stderr , "DEBUG 3, y is from %p to %p\n" , y , y + *n * *n * sizeof(double) ) ;
	free( y ) ; 
//fprintf( stderr , "DEBUG 3.5 qq is from %p to %p\n" , qq , qq + *n * *n * sizeof(double) ) ;
//fprintf( stderr , "DEBUG 4, idx ranges from %p to %p\n" , idx , idx + *n * sizeof(int) ) ; 
//fprintf( stderr , "DEBUG 4.5 " ) ; for( i = 0 ; i < *n + 1 ; i++ ) fprintf( stderr , "%i " , idx[i] ) ; fprintf( stderr , "\n" ) ; 
	free( idx ) ;   
//fprintf( stderr , "DEBUG 5\n" ) ;
	if( q != NULL ) 
		free( qq ) ; 
//fprintf( stderr , "DEBUG 6\n" ) ;
}

/*
void symmetricQrAlg ( double *mat , int *n , double *eps , double *eigVals , double *eigVecs , int *threads ) 
{
	// Convert to tri-diagonal form 
	double *tri = (double*) malloc( *n * *n * sizeof(double) ) ; 
	double *q = eigVecs ;
	symmToTri ( mat , n , tri , q , threads ) ; 
	double *a = eigVals ; 
	double *b = (double*) malloc( (*n - 1) * sizeof(double) ) ; 
	int i , j ; 
	for( i = 0 ; i < *n ; i++ ) 
	{
		a[i] = tri[ i + *n * i ] ; 
		if( i + 1 < *n ) 
			b[i] = tri[ i + 1 + *n * i ] ; 
	}
	
	// tri-diagonal QR Alg with Wilkinson shift 
	int m = *n ; 
	double d , c , s , x , y , err , t , prevT , w , z , tmp1 , tmp2 ; 
	double one = 1.0 ; 
	int DEBUGiter = 0 ; 
	while( m > 0 && DEBUGiter < 200 ) 
	{
DEBUGiter++ ; 
printf( "DEBUG, m: %i, iter: %i\n" , m , DEBUGiter ) ; 
		d = ( a[m-1] - a[m] ) / 2.0 ; // Wilkinson shift 
		if( d == 0.0 ) // Within tolerance ? 
			s = a[m] - fabs( b[m-1] ) ; 
		else 
			s = a[m] - b[m-1] * b[m-1] / ( d + copysign(one,d) * hypot( d , b[m-1] ) ) ; 
		x = a[0] - s ; 
		y = b[0] ; 
		for( i = 0 ; i < m - 1 ; i++ ) 
		{
			if( m > 1 ) 
				givens( x , y , &c , &s ) ; 
			else 
			{
				err = *eps + 1.0 ; 
				t = 1.111111111 ; 
				while( err > *eps ) 
				{
printf( "DEBUG, t: %f, err: %f\n" , t , err ) ; 
					prevT = t ; 
					t = t - ( (a[0]-a[1])*sin(2.0*t)/2.0 + b[0] * cos(2.0*t) ) / ( (a[0]-a[1])*cos(2.0*t) - 2.0*b[0]*sin(2.0*t) ) ; 
					err = fabs( prevT - t ) ; 
				}
			}
			w = c * x - s * y ; 
			d = a[i] - a[i+1] ; 
			z = ( 2.0 * c * b[i] + d * s ) * s ; 
			a[i] = a[i] - z ; 
			a[i+1] = a[i+1] + z ; 
			b[i] = d * c * s + (c*c - s*s) * b[i] ; 
			if( i > 0 ) 
				b[i-1] = w ; 
			if( i + 1 < m - 1 ) 
			{
				y = -s * b[i+1] ; 
				b[i+1] = c * b[i+1] ; 
			}
			for( j = 0 ; j < *n ; j++ ) 
			{
				tmp1 = c * q[ j , i ] - s * q[ j , i+1 ] ; 
				tmp2 = -s* q[ j , i ] + c * q[ j , i+1 ] ; 
				q[ j , i ] = tmp1 ; 
				q[ j , i+1 ] = tmp2 ; 
			}
		} 
printf( "DEBUG, err: %f\n" , fabs(b[m-1])/( fabs(a[m-1]) + fabs(a[m]) ) ) ; 
		if( fabs(b[m-1]) < *eps * ( fabs(a[m-1]) + fabs(a[m]) ) ) 
			m = m - 1 ; 
	}
	
	free( tri ) ; 
	free( b ) ; 
}
*/


























