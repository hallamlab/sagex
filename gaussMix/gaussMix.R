
dyn.load("gaussMix.so") 

c_mvdnorm = function( x , mu , sig ) 
{
	d = length(x) 
	out = .C( "mvdnorm" , as.double(x) , as.integer(d) , as.double(mu) , as.double(sig) , as.double(pi) , as.double(1.0) )[[6]] 
	return( exp(out) )  
}

# x : d X n matrix 
# p : k-length prob. vec. 
# mu : d X k matrix 
# sig : d X d X k k PSD matrices 
c_pGivenX = function( x , p , mu , sig ) 
{
	n = dim(x)[2] 
	d = dim(x)[1] 
	k = length(p) 
	out = matrix( 1:(k*n) , nrow=k , ncol=n ) 
	out = .C( "pGivenX" , as.double(x) , as.integer(n) , as.integer(d) , as.integer(k) , as.double(p) , as.double(mu) , as.double(sig) , as.double(pi) , as.double(out) )[[9]] 
	out = matrix( out , nrow=k , ncol=n ) 
	return( out ) 
}

c_eLogLik = function( x , p , mu , sig ) 
{
	n = dim(x)[2] 
	d = dim(x)[1] 
	k = length(p) 
	out = .C( "eLogLik" , as.double(x) , as.integer(n) , as.integer(d) , as.integer(k) , as.double(p) , as.double(mu) , as.double(sig) , as.double(pi) , as.double(1.0) )[[9]] 
	return( out ) 
}

# x : d X n matrix 
c_covMat = function( x )
{
	d = dim(x)[1] 
	n = dim(x)[2] 
	m = 1:d 
	cv = matrix( 1:(d*d) , nrow=d , ncol=d ) 
	tmp = .C( "covMat" , as.double(x) , as.integer(n) , as.integer(d) , as.double(m) , as.double(cv) ) 
	out = list() 
	out$mu = tmp[[4]] 
	out$sig = matrix( tmp[[5]] , nrow=d , ncol=d )  
	return( out ) 
}

c_boxMuller = function( n , m ) 
{
	out = matrix( 1:(n*m) , nrow=n , ncol=m ) 
	out = .C( "boxMuller" , as.integer(n) , as.integer(m) , as.double(pi) , as.double(out) )[[4]] 
	out = matrix( out , nrow=n , ncol=m ) 
	return( out ) 
}

# x : d X n matrix 
# k : integer > 0 
# sig output will be deformed 
c_init = function( x , k , threads=1 , eps=0.0000001 ) 
{
	n = dim(x)[2] 
	d = dim(x)[1] 
	p = 1:k 
	mu = matrix( 1:(d*k) , nrow=d , ncol=k ) 
	sig = matrix( 1:(d*d*k) , nrow=d , ncol=d*k ) 
	tmp = .C( "init" , as.double(x) , as.integer(n) , as.integer(d) , as.integer(k) , as.double(p) , as.double(eps) , as.integer(threads) , as.double(pi) , as.double(mu) , as.double(sig) ) 
	out = list() 
	out$p = tmp[[5]] 
	out$mu = matrix( tmp[[9]] , nrow=d , ncol=k ) 
	out$sig = matrix( tmp[[10]] , nrow=d , ncol=d*k ) 
	return( out ) 
}

c_nextP = function( pMat ) 
{
	k = dim(pMat)[1] 
	n = dim(pMat)[2] 
	out = 1:k 
	out = .C( "nextP" , as.double(pMat) , as.integer(n) , as.integer(1) , as.integer(k) , as.double(out) )[[5]] 
	return( out ) 
}

c_nextMu = function( x , pMat ) 
{
	d = dim(x)[1] 
	n = dim(x)[2] 
	k = dim(pMat)[1] 
	out = matrix( 1:(d*k) , nrow=d , ncol=k ) 
	out = .C( "nextMu" , as.double(x) , as.integer(n) , as.integer(d) , as.integer(k) , as.double(pMat) , as.double(out) )[[6]] 
	out = matrix( out , nrow=d , ncol=k ) 
	return( out ) 
}

c_nextSig = function( x , pMat , mu ) 
{
	d = dim(x)[1] 
	n = dim(x)[2] 
	k = dim(pMat)[1] 
	out = matrix( 1:(d*d*k) , nrow=d , ncol=d*k ) 
	out = .C( "nextSig" , as.double(x) , as.integer(n) , as.integer(d) , as.integer(k) , as.double(pMat) , as.double(mu) , as.double(out) )[[7]] 
	out = matrix( out , nrow=d , ncol=d*k ) 
	return( out ) 
}

c_fitMixture = function( x , k , eps=0.000001 , threads=1 , maxIter=100000 ) 
{
	n = dim(x)[2] 
	d = dim(x)[1] 
	p = 1:k 
	mu = matrix( 1:(d*k) , nrow=d , ncol=k ) 
	sig = matrix( 1:(d*d*k) , nrow=d , ncol=d*k ) 
	tmp = .C( "fitMixture" , as.double(x) , as.integer(n) , as.integer(d) , as.integer(k) , as.double(eps) , as.double(p) , as.double(mu) , as.double(sig) , as.integer(maxIter) , as.integer(threads) ) 
	out = list() 
	out$p = tmp[[6]] 
	out$mu = matrix( tmp[[7]] , nrow=d , ncol=k ) 
	out$sig = matrix( tmp[[8]] , nrow=d , ncol=d*k ) 
	return( out ) 
}

library(MASS) 

mvdnorm = function( x , mu , sig ) exp( -0.5*length(mu)*log(2*pi) - 0.5*log(det(sig)) - 0.5 * t(x - mu) %*% solve(sig) %*% (x - mu) )

# probability Choice i is j for each j 
# Choice i is denoted by choosing x == x[i,] 
pGivenX = function( p , x , mu , sig )
{
	den = function(j) ( p[j] * mvdnorm( x , mu[,j] , sig[[j]] ) ) 
	vals = sapply( 1:length(p) , den ) 
	return(  vals / sum(vals)  ) 
}

pMatGivenX = function( p , x , mu , sig ) 
{
	k = length(p) 
	n = dim(x)[1] 
	fun = function(i) pGivenX( p , t(x[i,]) , mu , sig ) 
	pVecs = lapply( 1:n , fun ) 
	out = as.matrix( pVecs[[1]] ) 
	if( n > 1 ) 
	{
		for( i in 2:n ) 
		{ out = cbind( out , pVecs[[i]] ) } 
	}
	return( out ) 
}

rDirichlet = function(alpha)
{
	n = length(alpha) 
	out = rgamma( n , alpha ) 
	return( out / sum(out) )  
}

eLogLik = function(x , p , mu , sig ) 
{
	k = length(p) 
	n = dim(x)[1] 
	out = 0.0 
	den = function(i,j) ( p[j] * mvdnorm( t(x[i,]) , mu[,j] , sig[[j]] ) ) 
	den3 = function(i)
	{
		den2 = function(j) den(i,j) 
		# return( sapply( 1:n , den2 ) ) # j is in 1:k, not 1:n !!! 
		return( sapply( 1:k , den2 ) ) 
	}
	cond = function(i) pGivenX( p , t(x[i,]) , mu , sig ) 
	eDen = function(i) sum( den3(i) * cond(i) ) 
	return( sum( sapply( 1:n , eDen ) ) ) 
}

nextMu = function( x , pMat ) 
{
	d = dim(x)[2] 
	n = dim(x)[1] 
	k = dim(pMat)[1] 
	# getAMuI = function(i,j) ( sum( t(x[,i]) * pMat[j,]) / sum( pMat[j,] ) ) 
	getAMu = function(j)
	{
		tmpFunc = function(i) ( sum( t(x[,i]) * pMat[j,]) / sum( pMat[j,] ) ) 
		return( sapply( 1:k , tmpFunc ) ) 
	}
	ll = lapply( 1:k , getAMu ) 
	out = ll[[1]] 
	if( k > 1 ) 
	{
		for( i in 2:k ) 
		{ out = cbind( out , ll[[i]] ) }
	} 
	return( out ) 
}

nextSig = function( x , mu1 , pMat ) 
{
	out = list()
	n = dim(x)[1] 
	d = dim(x)[2] 
	k = dim(mu1)[2] 
	for( j in 1:k ) 
	{
		const = sum( pMat[j,] ) 
		out[[j]] = matrix( rep(0,d*d) , nrow=d , ncol=d )  
		for( i in 1:n )
		{
			out[[j]] = out[[j]] + pMat[j,i] * ( t(x[i,]) - mu1[,j] ) %*% t( t(x[i,]) - mu1[,j] ) 
		}
		out[[j]] = out[[j]] / const 
	}
	return( out ) 
}

fitMix = function( x , k , eps=0.000001 , verbose=TRUE , maxIter=Inf ) 
{
	d = dim(x)[2] 
	n = dim(x)[1] 
	sig0 = cov(x) 
	mu0 = colSums(x)/n 
	tmp = eigen(sig0) 
	tmp = t(tmp$vectors) %*% diag( 0.5 * tmp$values ) %*% tmp$vectors 
	mu = t( mvrnorm( k , mu0 , tmp ) ) 
	sig = list() 
	for( i in 1:k ) 
	{ sig[[i]] = tmp } 
	# p = rDirichlet( rep( 1 , k ) ) 
	p = rep( 1/k , k ) 
	err = eps + 1.0 
	
	prev = eLogLik( x , p , mu , sig ) 
	count = 0 
	continue = TRUE  
	while( err > eps && count < maxIter && continue ) 
	{
		count = count + 1 
		
		pMat = pMatGivenX( p , x , mu , sig ) 
		p1 = rowSums( pMat ) / n 
		mu1 = nextMu( x , pMat ) 
		sig1 = nextSig( x , mu1 , pMat ) 
		
		post = eLogLik( x , p1 , mu1 , sig1 ) 
		if( post < prev ) # this should not happen 
		{ continue = FALSE }else{
			p = p1 
			mu = mu1 
			sig = sig1 
		}
		
		err = abs( prev - post ) 
		prev = post 
	}
	
	out = list() 
	out$p = p 
	out$mu = mu 
	out$sig = sig 
	out$iter = count 
	out$err = err 
	return( out ) 
}




