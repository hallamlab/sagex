
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <string.h>

// Appends two matrices 
// top : the matrix desired on the top 
// bot : the matrix desired on the bottom 
// rows1 : the number of rows in top 
// rows2 : the number of rows in bottom 
// cols : the number of columns in both matrices 
// out : enough space for (rows1 + rows2) * cols doubles 
void appendRows ( double *top , double *bot , int *rows1 , int *rows2 , int *cols , double *out ) ;  

// transpose a matrix 
// x : the matrix, also the output space 
// m : number of rows 
// n : number of columns 
// out : space for n X m doubles 
void transpose ( double *x , int *m , int *n , double *out ) ; 

// subtract matrix b from a 
// n : rows 
// m : cols 
// out : space for n X m doubles 
void subtract( double *a , double *b , int *n , int *m , double *out ) ;  

// Break a matrix into two matrices by dividing at a row 
// in : a column-major order matrix of n X m entries 
// div : number of rows in out1, from the upper part of 'in'  
// out1 , out2 : pre-allocated space 
void deappend ( double *in , int *n , int *m , int *div , double *out1 , double *out2 ) ; 

// Translates a matrix of ints to doubles 
// in : int matrix 
// n : rows 
// m : cols 
// out : a double matrix with rows X cols entries 
void intToDoubleMat ( int *in , int *n , int *m , double *out ) ; 

// Calculates the minor of a matrix 
// in : a column-indexed n X n matrix 
// n : pointer to an int storing the matrix dimension 
// i : the row to be removed 
// j : the column to be removed 
// out : an address for a single double to be stored 
void minorVal ( double *in , int *n , int *i , int *j , double *out ) ; 

// Calculates the determinant of matrix 
// in : a column-indexed n X n matrix 
// out : an address for a single double 
void det ( double *in , int *n , double *out ) ; 

// Calculates the adjoint matrix 
// in : a column-indexed n X n matrix 
// out : enough space for an n X n matrix 
void adj ( double *in , int *n , double *out ) ;  

/*
// Calculates the inverse of a matrix 
// in : a column-indexed n X n matrix 
// out : enough space for an n X n matrix 
void inv ( double *in , int *n , double *out ) ; 
*/

// calculates the inverse matrix of a positive-symmetric-definate matrix 
// in : a column-major order n X n PSD matrix to be inverted 
// out : space for n X n doubles 
// eps : convergence term 
// threads : number of available POSIX threads 
void invPsd( double *in , int *n , double *out , double *eps , int *threads ) ; 

// Calculates the product of matrices a & b 
// a : a column-indexed m X n matrix 
// b : a column-indexed n X p matrix 
// out : enough space for an m X p matrix 
void matProd ( double *a , double *b , int *m , int *n , int *p , double *out ) ; 

// Calculates the product of matrices a & b in parallel  
// a : a column-indexed m X n matrix 
// b : a column-indexed n X p matrix 
// out : enough space for an m X p matrix 
// threads : the number of threads to utilize 
void matProdP ( double *a , double *b , int *m , int *n , int *p , double *out , int *threads ) ; 

// Calculates the product of aT b
// a : a column-indexed n X m matrix 
// b : a column-indexed n X p matrix 
// out : enough space for an m X p matrix 
void matTrProd ( double *a , double *b , int *m , int *n , int *p , double *out ) ; 

// Calculates eigenvalues via the QR algorithm 
// in : a column-indexed n X n matrix 
// eps : convergence parameter, choose somthing small 
// val : space for n doubles, the eigen values 
// vec : space for an n X n matrix, the eigen vectors 
void qrEig ( double *in , int *n , double *eps , int *minIter , int *maxIter , double *val , double *vec ) ; 

// Calculates the QR decomposition 
// in : a column-indexed n X n matrix 
// q : space for an n X n matrix 
// r : space for an n X n matrix 
void qrDecomp ( double *in , int *n , double *q , double *r ) ; 

// Calculates an orthogonal basis via the Gram-Schmidt algorithm 
// in : a column-indexed n X n matrix 
// out : space for an n X n matrix 
void gramSchmidt ( double *in , int *n , double *out ) ; 

// Calculates the projection of vector x onto y 
void proj ( double *x , double *y , int *n , double *out ) ; 

// converts a symmetric matrix to a tri-diagonal matrix through a series of similar transforms 
// mat : a column-major n X n matrix 
// tri : enough space for n*n doubles 
// q : enough space for n*n doubles 
// void symmToTri ( double *mat , int *n , double *tri , double *q , int *threads ) ; 

// calculates eigen vecs & vals via power iteration on a symmetric matrix 
// accurracy may suffer slightly in exchange for more speed 
// mat : a symmetric n X n matrix 
// m : the number of desired eigen vectors 
// eps : tolerance, choose something small like e-14 
// val : output space, enough space for m eigen values 
// vec : output space, enough space for n X m doubles 
// threads : number of POSIX threads desired for the calculation 
void powerIteration( double *mat , int *n , int *m , double *eps , double *val , double *vec , int *threads ) ; 

// Calculates the Cholesky decomposition
// in : a column-major positive-definite n X n matrix 
// n : > 1 
// out : enough space for an n X n matrix 
void chol ( double *in , int *n , double *out ) ; 

