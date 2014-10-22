
dyn.load("crossValidate.so") 

crossValidate = function(k=1,p=85,a=2000)  
{
	gmN = 658928 + 1000
	sagN = 1000
	
	.C("crossValidate",as.integer(k),as.integer(p),as.integer(a))  
	try( {
	out = as.matrix(read.table("tmp6")) 
	.C( "cleanUp" ) 
	
	mat = matrix( 1:6 , nrow=2 , ncol=3 ) 
	rownames( mat ) = c( "is.sag" , "is.gm" ) 
	colnames( mat ) = c( "say.sag" , "say.gm" , "ppv" ) 
	
	sagObs = sum( substr(out[,1],2,5) == "Samp" ) 
	gmObs = dim(out)[1] - sagObs 
	
	mat[1,1] = sagObs 
	mat[2,1] = gmObs 
	mat[1,2] = sagN - sagObs 
	mat[2,2] = gmN - gmObs 
	
	mat[1,3] = mat[1,1] / sum( mat[,1] ) 
	mat[2,3] = mat[2,2] / sum( mat[,2] ) 
	
	outList = list() 
	outList$labels = out[,1] 
	outList$mat = mat
	
	return( outList ) 
	} )
	return( NULL ) 
}


