kVec = c( 0 , 5 , 15 , 16 , 17 , 18 , 19 , 20 , 25 , 30 ) 
correctHits1 = c( 1170 , 1170 , 1170 , 1169 , 1151 , 1056 , 931 , 884 , 687 , 671 ) 
hits1 = c( 1242 , 1242 , 1242 , 1238 , 1195 , 1072 , 934 , 885 , 687 , 673 ) 
correctHits2 = c( 1784 , 1784 , 1784 , 1784 , 1752 , 1581 , 1398 , 1311 , 1139 , 1103 ) 
hits2 = c( 2315 , 2314 , 2314 , 2301 , 2091 , 1719 , 1434 , 1326 , 1140 , 1106 ) 
falsePositives1 = hits1 - correctHits1 
falsePositives2 = hits2 - correctHits2 
trueNegatives1 = 658928 - falsePositives1 
trueNegatives2 = 658928 - falsePositives2 
falseNegatives1 = 2000 - correctHits1 
falseNegatives2 = 2000 - correctHits2 

plot( kVec , correctHits2 / hits2 , ylab="Percentage" , xlab="minimum shared subsequence length" , main="Positive Predictive Value (solid) and sensitivity (dashed) for EColi\nOne Cluster in Red, Two Clusters in Blue" , col="blue" , ylim = c(0.3,1.0) ) 
lines( kVec , correctHits2 / hits2 , type ='l' , lwd = 3 , col="blue" ) 
points( kVec , correctHits1 / hits1 , col="red" ) 
lines( kVec , correctHits1 / hits1 , type ='l' , lwd = 3 , col="red" )  
lines( kVec , correctHits2 / 2000 , lty=2 , type="o" , col = "blue" ) 
lines( kVec , correctHits1 / 2000 , lty=2 , type="o" , col = "red" )  

plot( kVec , trueNegatives1 / (trueNegatives1 + falseNegatives1) , ylab="Percentage" , xlab="minimum shared subsequence length" , main="Negative Predictive Value (solid) and specificity (dashed) for EColi\nOne Cluster in Red, Two Clusters in Blue" , col="red" , ylim=c( 0.9975,1 ) ) 
lines( kVec , trueNegatives1 / (trueNegatives1 + falseNegatives1) , lwd = 3 , col="red" ) 
points( kVec , trueNegatives2 / (trueNegatives2 + falseNegatives2) , col="blue" ) 
lines( kVec , trueNegatives2 / (trueNegatives2 + falseNegatives2) , lwd = 3 , col="blue" ) 
lines( kVec , trueNegatives1 / 658928 , lty=2 , type="o" , col="red" ) 
lines( kVec , trueNegatives2 / 658928 , lty=2 , type="o" , col="blue" ) 
