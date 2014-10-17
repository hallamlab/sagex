
dyn.load("kmeans.so") 

# x : d X n matrix 
c_kmeans = function(x,k, eps=0.00001 , maxIter=100000) 
{
	d = dim(x)[1] 
	n = dim(x)[2] 
	p = 1:k
	cov = matrix( 1:(d*d*k) , nrow=d , ncol=d*k ) 
	J = sample( 1:n , k ) 
	mu = x[,J] 
	tmp = .C( "kmeans" , as.double(x) , as.integer(d) , as.integer(n) , as.integer(k) , as.double(p) , as.double(mu) , as.double(cov) , as.double(eps) , as.integer(maxIter) ) 
	out = list() 
	out$p = tmp[[5]] 
	out$mu = matrix( tmp[[6]] , nrow=d , ncol=k )  
	out$cov = matrix( tmp[[7]] , nrow=d , ncol=d*k )
	return( out ) 
}

