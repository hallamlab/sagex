sag = read.table("ecoliSag.kmer")
gen = read.table("ecoliGen.kmer")
# plot metagenome 
plot( PC_2 ~ PC_0 , gen , pch=16 , col="grey" , cex=2 ) 
# plot Ecoli Genome 
points( PC_2 ~ PC_0 , gen[ gen$contigStatus == 2 ,] , pch=16 , col="red" , cex=2.1 ) 
# plot Ecoli SAG 
points( PC_2 ~ PC_0 , sag[ sag$contigStatus == 2 ,] , pch=16 , col="darkred" , cex=2.1 )  
