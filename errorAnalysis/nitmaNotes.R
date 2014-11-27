
kVec = c( 0 , 5 , 15 , 16 , 17 , 18 , 19 , 20 , 25 , 30 , 40 , 50 , 60 , 70 , 80 ,90 , 100 ) 
trueHits1 = c( 1085 , 1085 , 1084 , 1088 , 1084 , 1073 , 1060 , 1054 , 1243 , 1258 , 1279 , 1301 , 1312 , 1319 , 1319 , 1321 , 1319 ) 
trueHits2 = c( 1937 , 1937 , 1937 , 1935 , 1934 , 1932 , 1928 , 1919 , 1922 , 1935 , 1918 , 1916 , 1923 , 1926 , 1927 , 1925 , 1925 ) 
hits1 = c( 4755 , 4750 , 4749 , 4745 , 4389 , 3380 , 2293 , 1427 , 1598 , 1516 , 1425 , 1363 , 1349 , 1339 , 1328 , 1324 , 1320 ) 
hits2 = c( 23711 , 23695 , 23694 , 23518 , 21127 , 15017 , 8953 , 5303 , 3234 , 2956 , 2429 , 2150 , 2046 , 1982 , 1951 , 1932 , 1926 ) 
falsePositives1 = hits1 - trueHits1 
falsePositives2 = hits2 - trueHits2 
trueNegatives1 = 658928 - falsePositives1 
trueNegatives2 = 658928 - falsePositives2 
falseNegatives1 = 2000 - trueHits1 
falseNegatives2 = 2000 - trueHits2 


plot( kVec , trueHits1 / hits1 , ylab="Percentage" , xlab="minimum shared subsequence length" , main="Positive Predictive Value (solid) and sensitivity (dashed) of idealized Nitma\nOne cluster in red, two clusters in blue" , col="red" , ylim=c(0,1) ) 
lines( kVec , trueHits1 / hits1 , lwd=3 , col="red" ) 
points( kVec , trueHits2 / hits2 , col="blue" ) 
lines( kVec , trueHits2 / hits2 , lwd=3 , col="blue" ) 
lines( kVec , trueHits1 / 2000 , lty=2 , type="o" , col="red" ) 
lines( kVec , trueHits2 / 2000 , lty=2 , type="o" , col="blue" ) 

plot( kVec , trueNegatives1/(trueNegatives1 + falseNegatives1) , ylab="Percentage" , xlab="minimum shared subsequence length" , main="Negative Predictive Value (solid) and Sensitivity (dashed) of idealized Nitma\nOne cluster in red, two clusters in blue" , col="red" , ylim=c(0.95,1) ) 
lines( kVec , trueNegatives1/(trueNegatives1 + falseNegatives1) , lwd=3 , col="red" ) 
points( kVec , trueNegatives2/(trueNegatives2 + falseNegatives2) , col="blue" ) 
lines( kVec ,trueNegatives2/(trueNegatives2 + falseNegatives2) , lwd=3 , col="blue" ) 
lines( kVec , trueNegatives1/658928 , lty=2 , type="o" , col="red" ) 
lines( kVec , trueNegatives2/658928 , lty=2 , type="o" , col="red" ) 

