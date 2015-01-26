#!/bin/bash

#CC = g++
#CFLAGS = -Wall
#DEPS = parse_Blastoutput.h
#OBJ = SAG_extrapolator.o, parse_Blastoutput.o
#
#%.o: %.cpp $(DEPS)
#	$(CC) $(CFLAGS) -c -o $@ $<
#
#SAG_extrapolator: $(OBJ)
#	g++ $(CFLAGS) -o $@ $^

cp ./posixMat/* ./build/
cp ./chopKmers/* ./build/
cp ./gammaDist/* ./build/
cp ./classify/* ./build/
cp ./interface/* ./build/
cp ./gaussMix/* ./build/
cp ./gmmPval/* ./build/
cp ./kmeans/* ./build/ 
cp ./kbag/* ./build/ 

rm ./build/*.o
rm ./build/*.so
rm ./build/*.R

g++ -g -o sagex ./build/SAG_extrapolator.cpp ./build/parse_Blastoutput.cpp ./build/FastaParser.cpp ./build/Genome.cpp ./build/metBagger.cpp ./build/count.c ./build/kmers.c ./build/minLength.c ./build/recordTasks.c ./build/gammaDist.c ./build/matrix.c ./build/stat.c ./build/gmmPval.c ./build/gaussMix.c ./build/kmeans.c ./build/classify.c ./build/kbag.c -pthread


