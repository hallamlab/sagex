
dyn.load("stat.so") 

colStan = function( x ) 
{
	n = dim(x)[1] 
 	m = dim(x)[2] 
	tmp = matrix( 1:(n*m) , nrow=n , ncol=m ) 
	tmp = .C( "colStandardize" , as.double(x) , as.integer(n) , as.integer(m) , as.double(tmp) )[[4]] 
	tmp = matrix( tmp , nrow=n , ncol=m ) 
	return( tmp ) 
}

corrMat = function( x ) 
{
	n = dim(x)[1] 
	m = dim(x)[2] 
	tmp = matrix( 1:(m*m) , nrow=m , ncol=m ) 
	tmp = .C( "corrMat" , as.double(x) , as.integer(n) , as.integer(m) , as.double(tmp) )[[4]] 
	tmp = matrix( tmp , nrow=m , ncol=m ) 
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

