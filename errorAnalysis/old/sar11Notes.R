
kVec = c( 0 , 5 , 15 , 16 , 17 , 18 , 19 , 20 ,25 , 30 , 40 , 50 , 60 , 70 , 80 , 90 , 100 ) 
trueHits1 = c( 787 , 786 , 787 , 790 , 780 , 750 , 615 , 479 , 325 , 244 , 114 , 54 , 32 ,32 , 28 , 16 , 14 ) 
trueHits2 = c( 1893 , 1893 , 1893 , 1893 , 1875 , 1736 , 1400 , 1040 , 583 , 421 , 203 , 107 , 73 , 67 , 53 , 36 , 36 ) 
hits1 = c( 4355 , 4350 , 4356 , 4364 , 4178 , 3383 , 2237 , 1372 , 633 , 614 , 352 , 223 , 152 , 152 , 157 , 118 , 99 ) 
hits2 = c( 16113 , 16102 , 16105 , 16094 , 15426 , 12608 , 8350 , 4999 , 1557 , 1228 , 762 , 539 , 415 , 356 , 288 , 217 , 190 ) 
falsePositives1 = hits1 - trueHits1 
falsePositives2 = hits2 - trueHits2 
trueNegatives1 = 658928 - falsePositives1 
trueNegatives2 = 658928 - falsePositives2 
falseNegatives1 = 2000 - trueHits1 
falseNegatives2 = 2000 - trueHits2 

plot( kVec , trueHits1 / hits1 , ylab="Percentage" , xlab="minimum shared subsequence length" , main="Positive Predictive Value (solid) and Sensitivity (dashed) for Pelagibacter\nOne cluster in red, two clusters in blue" , col="red" , ylim=c(0,1) ) 
lines( kVec , trueHits1 / hits1 , lwd=3 , col="red" ) 
points( kVec , trueHits2 / hits2 , col="blue" ) 
lines( kVec , trueHits2 / hits2 , lwd=3 , col="blue" ) 
lines( kVec , trueHits1 / 2000 , lty=3 , type="o" , col="red" ) 
lines( kVec , trueHits2 / 2000  , lty=3 , type="o" , col="blue" ) 

plot( kVec , trueNegatives1 / (trueNegatives1 + falseNegatives1) , ylab="Percentage" , xlab="minimum shared subsequence length" , main="Negative Predictive Value (solid) and Specificity (dashed) for Pelagibacter\nOne cluster in red, two clusters in blue" , col="red" , ylim=c(0.975,1) ) 
lines( kVec , trueNegatives1 / (trueNegatives1 + falseNegatives1) , lwd=3 , col="red" ) 
points( kVec , trueNegatives2 / (trueNegatives2 + falseNegatives2) , col="blue" ) 
lines( kVec , trueNegatives2 / (trueNegatives2 + falseNegatives2) , lwd=3 , col="blue" ) 
lines( kVec , trueNegatives1 / 658928 , lty=3 , type="o" , col="red" ) 
lines( kVec , trueNegatives2 / 658928 , lty=3 , type="o" , col="blue" ) 

