
dyn.load("gmmPval.so") 

qBootGMM = function( q , p , mu , sig , n=10000 , eps=0.00001 , threads=1 ) 
{
	d = dim(mu)[1] 
	k = dim(mu)[2] 
	sigg = sig[[1]] 
	if( k > 1 ) 
	{
		for( i in 2:k ) 
		{ sigg = cbind( sigg , sig[[i]] ) } 
	}
	out = .C( "getGMMQuantile" , as.double(q) , as.double(p) , as.double(mu) , as.double(sigg) , as.integer(n) , as.integer(d) , as.integer(k) , as.double(eps) , as.integer(threads) , as.double(1.0) )[[10]] 
	return( out ) 
}

# simulates data from a GMM of k gaussians 
# n : number of simulants 
# mu : a k X d matrix of mu vectors 
# p : a k-length array probability vector 
# sig : a list of k d X d PSD matrices 
simulateGMM = function( n , mu , p , sig , eps=0.00001 , threads=1 ) 
{
	k = dim(mu)[1] 
	d = dim(mu)[2] 
	sigg = sig[[1]] 
	if( k > 1 ) 
	{
		for( i in 2:k ) 
		{ sigg = cbind( sigg , sig[[i]] ) } 
	}
	out = matrix( 1:(d*n) , nrow=d , ncol=n ) 
	out = .C( "simulateGMM" , as.integer(n) , as.integer(d) , as.integer(k) , as.double(p) , as.double(mu) , as.double(sigg) , as.double(eps) , as.integer(threads) , as.double(out) )[[9]] 
	out = matrix( out , nrow=d , ncol=n ) 
	return( out ) 
}

matMuller = function( n , m ) 
{
	out = matrix( 1:(n*m) , nrow=n , ncol=m ) 
	out = .C( "matMuller" , as.integer(n) , as.integer(m) , as.double(out) )[[3]] 
	out = matrix( out , nrow=n , ncol=m ) 
	return( out ) 
}

#quickSort = function( x ) 
#{
#	n = length(x) 
#	out = .C( "quickSort" , as.double(x) , as.integer(n) )[[1]]  
#	return( out ) 
#}

