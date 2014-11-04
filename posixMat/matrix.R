
dyn.load("matrix.so") 

a = matrix( rnorm(10) , nrow=5, ncol=2 ) 
b = matrix( rnorm(8) , nrow=2 , ncol=4 ) 
c = matrix( rpois(25,25) , nrow=5 , ncol=5 ) 
cc = cov(c)
g = matrix( rnorm(100) , nrow=20 , ncol = 5 ) 
gg = cov(g) 

householder = function(x,k) 
{
	n = dim(x)[1] 
	out = c( rep(0,k) , x[(k+1):n,k] ) 
	ttl = sqrt( sum( out * out ) ) 
	out[k+1] = out[k+1] + sign(out[k+1]) * ttl 
	ttl = sqrt( sum( out * out ) ) 
	out = out / ttl 
	u = out 
	v = as.vector( 2 * x %*% u - 2 * u * ( t(u) %*% x %*% u ) )  
	out = list() 
	out$u = u 
	out$v = v 
	return( out ) 
} 

matHouseholder = function(x,k) 
{
	n = dim(x)[1] 
	tmp = .C( "householder" , as.double(x) , as.integer(n) , as.integer(k) , as.double(1:n) , as.double(1:n) ) 
	u = tmp[[4]] 
	v = tmp[[5]] 
	out = list() 
	out$u = u 
	out$v = v 
	return( out )  
}

quickSort = function( x ) 
{
        n = length(x) 
        tmp = .C( "quickSort" , as.double(x) , as.integer(n) , as.integer(1:n) ) 
        out = list() 
	out$x = tmp[[1]] 
	out$idx = tmp[[3]] 
	return( out ) 
}

givens = function(x,y) 
{
	r = sqrt( x*x + y*y ) 
	c = x/r 
	s = -y/r 
	out = list() 
	out$c = c 
	out$s = s 
	return( out )  
}

matGivens = function(x,y) 
{
	tmp = .C( "givens" , as.double(x) , as.double(y) , as.double(1.0) , as.double(1.0) ) ; 
	out = list() 
	out$c = tmp[[3]] 
	out$s = tmp[[4]] 
	return( out )  
}

matPsdEig = function( x , eps=0.0000000000000001 ) 
{
	n = dim(x)[1] 
	q = matrix( 1:(n*n) , nrow=n , ncol=n ) 
	v = 1:n
	tmp = .C( "psdEig" , as.double(x) , as.integer(n) , as.double(eps) , as.double(q) , as.double(v) ) 
	out = list() 
	out$q = matrix( tmp[[4]] , nrow=n , ncol=n )
	out$v = tmp[[5]] 
	return( out ) 
}

psdEig = function( x , eps=0.00000000000000001 )  
{
	n = dim(x)[1]  
	y = x 
	Q = diag( rep(1,n) ) 
	# convert to tridiag 
	for( i in 1:(n-1) ) 
	{
		hh = householder(y,i) 
		Q = Q %*% ( diag(1,n) - 2.0 * hh$u %*% t(hh$u) ) 
		y = y - hh$u %*% t(hh$v) - hh$v %*% t(hh$u) 
	}
	
	for( m in n:2 ) 
	{
		err = eps + 1 
		iter = 0 
		while( err > eps ) 
		{
			iter = iter + 1 
			# Wilkinson shift 
			d = ( y[m-1,m-1] - y[m,m] ) / 2 
			s = NULL 
			if( d == 0 ) { s = y[m,m] - abs(y[m-1,m]) } else { s = y[m,m] - (y[m-1,m]^2) / ( d + sign(d) * sqrt( d*d + (y[m-1,m]^2) ) ) } 
			xx = y[1,1] - s 
			yy = y[1,2] 
			for( k in 1:(m-1) ) 
			{
				c = NULL 
				s = NULL 
				if( m > 2 ) 
				{
					g = givens(xx,yy) 
					c = g$c 
					s = g$s 
				}
				else
				{
					t = 0.5 * atan( 2 * y[1,2] / (y[2,2] - y[1,1]) ) 
					c = cos(t) 
					s = sin(t) 
				} 
				w = c*xx - s*yy 
				d = y[k,k] - y[k+1,k+1] 
				z = s * ( 2 * c * y[k,k+1] + d * s ) 
				y[k,k] = y[k,k] - z 
				y[k+1,k+1] = y[k+1,k+1] + z 
				y[k,k+1] = d * c * s + (c*c - s*s) * y[k,k+1] 
				y[k+1,k] = y[k,k+1] 
				xx = y[k,k+1] 
				if( k > 1 ) { y[k-1,k] = w ; y[k,k-1] = w } 
				if( k < m - 1 ){ yy = -s * y[k+1,k+2] ; y[k+1,k+2] = c * y[k+1,k+2] ; y[k+2,k+1] = y[k+1,k+2] } 
				Q[,k:(k+1)] = Q[,k:(k+1)] %*% cbind( c(c,-s) , c(s,c) ) 
				# err = abs(y[m-1,m]) / ( abs( y[m-1,m-1] ) + abs( y[m,m] ) ) 
			}
			err = abs(y[m-1,m]) / ( abs( y[m-1,m-1] ) + abs( y[m,m] ) ) 
		}
	}
	
	out = list() 
	out$v = diag(y) 
	idx = sort( out$v , index.return=T , decreasing=T )$ix 
	out$v = out$v[ idx ] 
	out$q = Q[ , idx ] 
	return( out )  
} 

library(MASS)
bigMat = function( n )
{
	cc = rWishart( 1 , n , diag( rep(1,n) ) )[,,1]  
	return( cov( mvrnorm( 3*n , rep(0,n) , cc ) ) ) 
}

matTr = function( x ) 
{
	m = dim(x)[1] 
	n = dim(x)[2] 
	out = matrix( 1:(m*n) , nrow=n, ncol=m) 
	out = .C( "transpose" , as.double(x) , as.integer(m) , as.integer(n) , as.double(out) )[[4]] 
	out = matrix( out , nrow = n , ncol = m ) 
	return( out ) 
}

matProd = function( a , b ) 
{
	m = dim(a)[1] 
	n = dim(a)[2] 
	if( n != dim(b)[1] )
	{ return(NULL) } 
	p = dim(b)[2] 
	tmp = matrix( 1:(m*p) , nrow=m , ncol=p ) 
	tmp = .C( "matProd" , as.double(a) , as.double(b) , as.integer(m) , as.integer(n) , as.integer(p) , as.double(tmp) )[[6]]  
	tmp = matrix( tmp , nrow=m , ncol=p ) 
	return( tmp ) 
}

matProdP = function( a , b , threads=1 )
{
        m = dim(a)[1]
        n = dim(a)[2]
        if( n != dim(b)[1] || threads < 1 )
        { return(NULL) }
        p = dim(b)[2]
        tmp = matrix( 1:(m*p) , nrow=m , ncol=p )
        tmp = .C( "matProdP" , as.double(a) , as.double(b) , as.integer(m) , as.integer(n) , as.integer(p) , as.double(tmp) , as.integer(threads) )[[6]]
        tmp = matrix( tmp , nrow=m , ncol=p )
        return( tmp )
}

matTrProd = function( a , b ) 
{
	m = dim(a)[1] 
	if( m != dim(b)[1] ) 
	{ return( NULL ) } 
	n = dim(a)[2] 
	p = dim(b)[2] 
	tmp = matrix( 1:(n*p) , nrow=n , ncol=p ) 
	tmp = .C( "matTrProd" , as.double(a) , as.double(b) , as.integer(m) , as.integer(n) , as.integer(p) , as.double(tmp) )[[6]] 
	tmp = matrix( tmp , nrow = m , ncol = p ) 
	return( tmp ) 
}

matCov = function( x )
{
	n = dim(x)[1] 
	m = dim(x)[2] 
	mu = 1:m 
	sig = matrix( 1:(m*m) , nrow=m , ncol=m ) 
	tmp = .C( "covMat" , as.double(x) , as.integer(n) , as.integer(m) , as.double(mu) , as.double(sig) ) 
	mu = tmp[[4]] 
	sig = matrix( tmp[[5]] , nrow=m , ncol=m )
	out = list()
	out$mu = mu 
	out$sig = sig 
	return( out ) 
}

matApp = function( a , b ) 
{
	n = dim(a)[1] 
	m = dim(b)[1] 
	cols = dim(a)[2] 
	if( cols != dim(b)[2] ) 
	{ return(NULL) } 
	tmp = matrix( 1:((n+m)*cols) , nrow = n + m , ncol = cols ) 
	tmp = .C( "appendRows" , as.double(a) , as.double(b) , as.integer(n) , as.integer(m) , as.integer(cols) , as.double(tmp) )[[6]] 
	tmp = matrix( tmp , nrow= n + m , ncol= cols ) 
	return( tmp )  
}

matDeapp = function( a , div ) 
{
	n = dim(a)[1] 
	m = dim(a)[2] 
	tmp1 = matrix( 1:(div*m) , nrow=div , ncol = m ) 
	tmp2 = matrix( 1:((n-div)*m) , nrow=n-div , ncol = m ) 
	tmp = .C( "deappend" , as.double(a) , as.integer(n) , as.integer(m) , as.integer(div) , as.double(tmp1) , as.double(tmp2) ) 
	tmp1 = matrix( tmp[[5]] , nrow = div , ncol = m ) 
	tmp2 = matrix( tmp[[6]] , nrow = n - div , ncol = m ) 
	out = list()
	out$hi = tmp1 
	out$low = tmp2 
	return( out ) 
}

matDet = function( a ) 
{
	n = dim(a)[1]
	if( n != dim(a)[2] ) 
	{ return( NULL ) } 
	out = .C( "det" , as.double(a) , as.integer(n) , as.double(5) )[[3]] 
	return( out ) 
}

matInvPsd = function( a , eps=0.00000001 , threads=1 ) 
{
	n = dim(a)[1] 
	if( n != dim(a)[2] ) 
	{ return( NULL ) } 
	if( matDet(a) < tol )
	{
		print( "Matrix computationally singular" ) 
		return( matDet(a) ) 
	}
	tmp = matrix( 1:(n*n) , nrow = n , ncol = n ) 
	tmp = .C( "invPsd" , as.double(a) , as.integer(n) , as.double(tmp) , as.double(eps) , as.integer(threads) )[[3]] 
	tmp = matrix( tmp , nrow = n , ncol = n ) 
	return( tmp ) 
}

matQr = function( a , tol=0.00000000001 ) 
{
	n = dim(a)[1] 
	if( n != dim(a)[2] ) 
	{ return( NULL ) } 
	if( matDet(a) < tol ) 
	{
		print( "Matrix computationally singular" ) ; 
		return( matDet(a) ) 
	}
	q = matrix( 1:(n*n) , nrow=n , ncol=n ) 
	r = q
	tmp = .C( "qrDecomp" , as.double(a) , as.integer(n) , as.double(q) , as.double(r) )  
	q = tmp[[3]] 
	r = tmp[[4]] 
	out = list()
	out$q = matrix(q,nrow=n,ncol=n)  
	out$r = matrix(r,nrow=n,ncol=n)  
	return( out ) 
}

matEig = function( a , tol=0.00000000001 , maxIter=1000 , minIter=10 ) 
{
	n = dim(a)[1] 
	if( n != dim(a)[2] ) 
	{ return(NULL) } 
	# TODO speed up matDet before running this 
	#if( matDet(a) < tol )
	#{
	#	print( "Matrix computationally singular" ) 
	#	return( matDet(a) ) 
	#}
	vec = matrix( 1:(n*n) , nrow=n , ncol=n ) 
	val = 1:n
	tmp = .C( "qrEig" , as.double(a) , as.integer(n) , as.double(tol) , as.integer(minIter) , as.integer(maxIter) , as.double(val) , as.double(vec) ) 
	val = tmp[[6]] 
	vec = matrix( tmp[[7]] , nrow=n , ncol=n ) 
	out = list() ; 
	out$val = val 
	out$vec = vec 
	out$iter = tmp[[5]] 
	return( out ) 
}

matProj = function( x , y ) # project x onto y 
{
	n = length(x) 
	if( n != length(y) || n == 0 ) 
	{ return( NULL ) } 
	out = .C( "proj" , as.double(x) , as.double(y) , as.integer(n) , as.double(1:n) )[[4]] 
	return( out ) 
}

matGs = function( x ) 
{
	n = dim(x)[1] 
	if( n != dim(x)[2] ) 
	{ return(NULL) } 
	tmp = matrix( 1:(n*n) , nrow=n , ncol=n ) 
	tmp = .C( "gramSchmidt" , as.double(x) , as.integer(n) , as.double(tmp) )[[3]]  
	tmp = matrix( tmp , nrow=n , ncol=n ) 
	return( tmp ) 
}

matChol = function( x ) 
{
	n = dim(x)[1] 
	if( n != dim(x)[2] )
	{ return(NULL) } 
	if( sum( x == t(x) ) != n*n ) 
	{ return(NULL) } 
	tmp = matrix( 1:(n*n) , nrow=n , ncol=n ) 
	tmp = .C( "chol" , as.double(x) , as.integer(n) , as.double(tmp) )[[3]] 
	tmp = matrix( tmp , nrow=n , ncol=n ) 
	return( tmp ) 
}

matInvPsd = function( x , eps=0.000000001 , threads=1 ) 
{
	n = dim(x)[1] 
	if( n != dim(x)[2] ) 
	{ return(NULL) } 
	tmp = matrix( 1:(n*n) , nrow=n , ncol=n ) 
	tmp = .C( "invPsd" , as.double(x) , as.integer(n) , as.double(tmp) , as.double(eps) , as.integer(threads) )[[3]]
	tmp = matrix( tmp , nrow=n , ncol=n ) 
	return( tmp )  
}

matSymmToTri = function( x , threads=1 ) 
{
	n = dim(x)[1] 
	if( n != dim(x)[2] || threads < 1 ) 
	{ return(NULL) } 
	tri = matrix( 1:(n*n) , nrow = n , ncol = n ) 
	q = matrix( 1:(n*n) , nrow=n , ncol=n ) 
	tmp = .C( "symmToTri" , as.double(x) , as.integer(n) , as.double(tri) , as.double(q) , as.integer(threads) ) 
	tri = matrix( tmp[[3]] , nrow=n , ncol=n ) 
	q = matrix( tmp[[4]] , nrow= n , ncol=n ) 
	out = list() 
	out$tri = tri 
	out$q = q 
	return( out ) 
}

matSymmQR = function( x , eps=0.00001 , threads=1 ) 
{
	n = dim(x)[1] 
	if( n != dim(x)[2] ) 
	{ return( NULL ) } 
	vecs = matrix( 1:(n*n) , nrow=n , ncol=n ) 
	vals = 1:n 
	tmp = .C( "symmetricQrAlg" , as.double(x) , as.integer(n) , as.double(eps) , as.double(vals) , as.double(vecs) , as.integer(threads) ) ; 
	vecs = matrix( tmp[[5]] , nrow=n , ncol=n ) 
	vals = tmp[[4]]   
	out = list() 
	out$vals = vals 
	out$vecs = vecs 
	return( out ) 
}

matPwrTri = function( tri , x ) 
{
	n = dim(tri)[1] 
	tmp = .C( "powIterTridiag" , as.double(tri) , as.integer(n) , as.double(x) , as.double(1:n) , as.double(1) )
	vec = tmp[[4]] 
	val = tmp[[5]] 
	out = list() 
	out$vec = vec 
	out$val = val 
	return( out ) 
}

matRmDim = function( x , y ) # remove y from x 
{
	n = dim(x)[1] 
	out = matrix( 1:(n*n) , nrow=n , ncol=n) 
	out = .C( "removeDimension" , as.double(x) , as.integer(n) , as.double(y) , as.double(out) )[[4]] 
	out = matrix( out , nrow=n , ncol=n ) 
	return( out ) 
}

matGetEigenTri = function( tri , eps=0.00000000000001 )
{
	n = dim(tri)[1] 
	tmp = .C( "getEigenPowerItr" , as.double(tri) , as.integer(n) , as.double(eps) , as.double(1) , as.double(1:n) )  
	out = list() 
	out$val = tmp[[4]] 
	out$vec = tmp[[5]] 
	return( out ) 
}

matPwrSym = function( mat , m , eps=0.00000000000001 , threads=1 ) 
{
	n = dim(mat)[1] 
	tmp = .C( "powerIteration" , as.double(mat) , as.integer(n) , as.integer(m) , as.double(eps) , as.double(1:m) , as.double(1:(n*m)) , as.integer(threads) )  
	out = list() 
	out$vals = tmp[[5]] 
	out$vec = matrix( tmp[[6]] , nrow=n , ncol=m ) 
	return( out ) 
}







