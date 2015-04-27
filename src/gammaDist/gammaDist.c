
#include <math.h>
#include <stdio.h>

void density ( double *beta , double *alpha , double *x , double *out )
{
	*out = exp( *alpha * log( *beta ) + (*alpha-1.0) * log( *x ) - *beta * *x - lgamma(*alpha)  ) ;  
}


// Integrates the gamma density from 0 to x via Riemann sum with eps-sized bins
// up : > 0 if integrating from + side of bins, otherwise - 
void intGamma ( double *beta , double *alpha , double *x , double *eps , int *up , double *out )
{
	*out = 0.0 ; 
	double y , tmp , pos ; 
	for( y = 2.0 * *eps ; y < *x ; y += *eps ) 
	{
		if( *up > 0 )
			pos = y ; 
		else
			pos = y - *eps ; 
		density( beta , alpha , &pos , &tmp ) ; 
		*out += tmp * *eps ; 
	}
}

void deriv ( double *beta , double *alpha , double *x , double *out ) 
{
	density( beta , alpha , x , out ) ; 
	*out = *out * ( (*alpha - 1.0) / (*x) - *beta ) ; 
}

// Invert the gamma distribution with Newton-Raphson 
// alpha : distn param 1 
// beta : distn param 2 
// p : desired cumulative probability 
void nrGammaInv ( double *beta , double *alpha  , double *p , double *eps , double *init , int *maxIter , double *out ) 
{
	*out = *alpha / *beta ; // E[X] 
	if( init != NULL ) 
		*out = *init ;  
	int i ; 
	double err = *eps + 1.0 ; 
	double prevVal , num , den ; 
	int zero = 0 ; 
	for( i = 0 ; err > *eps && i < *maxIter ; i++ )	
	{
printf( "DEBUG: %f\n" , *out ) ; 
		prevVal = *out ; 
		intGamma( beta , alpha , out , eps , &zero , &num ) ; 
		num = num - *p ; 
		deriv( beta , alpha , out , &den ) ; 
		*out = *out - num / den ; 
		err = fabs( *out - prevVal ) ; 
	}
}


// Dyadic inversion for gamma distirbution 
// alpha : distn param 1 
// beta : distn param 2 
// p : desired cumulative probability 
void dyadicInvGam ( double *beta , double *alpha , double *p , double *eps , int *maxIter , double *out ) 
{
	*out = *alpha / *beta ; // E[X] 
	double prevVal ; 
	double err = *eps + 1.0 ; 
	int zero = 0 ; 
	double tmp ; 
	int flag ; 
	intGamma ( beta , alpha , out , eps , &zero , &tmp ) ; 
	if( tmp < *p )
	{
		flag = 1 ; // out must increase 
		*out = *out * 2.0 ; 
	}
	else
	{
		flag = -1 ; // out must decrease 
		*out = *out / 2.0 ; 
	}
	double upper , lower ; 
	int i ; 
	for( i = 0 ; err > *eps && i < *maxIter ; i++ ) 
	{
		prevVal = *out ; 
		intGamma ( beta , alpha , out , eps , &zero , &tmp ) ; 
		if( flag > 0 ) // searching up 
		{
			if( tmp > *p ) // upper found 
			{
				upper = *out ; 
				lower = *out / 2.0 ; 
				*out = (upper + lower) / 2.0 ; 
				flag = 0 ; 
			}
			else
				*out = *out * 2.0 ; 
		}
		if( flag < 0 ) // searching down 
		{
			if( tmp < *p ) // lower found 
			{
				lower = *out ; 
				upper = *out * 2.0 ; 
				*out = (upper + lower) / 2.0 ; 
				flag = 0 ; 
			}
			else
				*out = *out / 2.0 ; 
		}
		if( flag == 0 ) // searching dyadically 
		{
			if( tmp < *p ) // search up 
			{
				lower = *out ; 
				*out = (upper + lower) / 2.0 ; 
			}
			else // search down 
			{
				upper = *out ; 
				*out = (upper + lower) / 2.0 ; 
			}
			err = fabs( *out - prevVal ) ; 
		}
	}
}























