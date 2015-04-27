
dyn.load("gammaDist.so") 

dist = function( beta , alpha , x , eps=0.001 ) 
{
	zero = 0 
	out = .C( "intGamma" , as.double(beta) , as.double(alpha) , as.double(x) , as.double(eps) , as.integer(zero) , 1 )[[6]] 
	return( out ) 
}

invGam = function( beta , alpha , p , eps=0.0001 , maxIter=1000 ) 
{
	out = .C( "dyadicInvGam" , as.double(beta) , as.double(alpha) , as.double(p) , as.double(eps) , as.integer(maxIter) , 1 )[[6]] 
	return( out ) 
}

