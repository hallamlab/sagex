
#===========================
# SAGex Makefile
# Connor Morgan-Lang 2015
#===========================

SHELL := /bin/bash -e

export OBJ_DIR = obj#For the compiled object files
export SRC_DIR = src#Where the c++ files live

CC=g++
CFLAGS=-Wall -c
INC=-I$(SRC_DIR)/classify -I$(SRC_DIR)/posixMat -I$(SRC_DIR)/gammaDist -I$(SRC_DIR)/chopKmers -I$(SRC_DIR)/gaussMix -I$(SRC_DIR)/gmmPval -I$(SRC_DIR)/kbag -I$(SRC_DIR)/kmeans -I$(SRC_DIR)/interface

OBJECTS = $(OBJ_DIR)/posixMat.o \
	$(OBJ_DIR)/kmeans.o \
	$(OBJ_DIR)/kbag.o \
	$(OBJ_DIR)/gmmPval.o \
	$(OBJ_DIR)/gaussMix.o \
	$(OBJ_DIR)/gammaDist.o \
	$(OBJ_DIR)/classify.o \
	$(OBJ_DIR)/count.o \
	$(OBJ_DIR)/kmers.o \
	$(OBJ_DIR)/recordTasks.o \
	$(OBJ_DIR)/minLength.o \
	$(OBJ_DIR)/FastaParser.o \
	$(OBJ_DIR)/interface.o 

all: Notify_beginning $(OBJ_DIR) sagex

Notify_beginning:
	@echo "****************************"
	@echo "Making SAGex!"
	@echo "****************************"

$(OBJ_DIR):
	@mkdir -p $@

#sagex: $(OBJECTS)
#	$(CC) -o sagex $(OBJECTS)

#A utils directory would be handy in the future to make the sub-src directories tidier

$(OBJ_DIR)/posixMat.o: $(SRC_DIR)/posixMat/matrix.c 
	$(CC) $(CFLAGS) src/posixMat/matrix.c src/posixMat/matrix.h -pthread -o $@

$(OBJ_DIR)/kmeans.o: $(SRC_DIR)/kmeans/kmeans.c 
	$(CC) $(CFLAGS) $(SRC_DIR)/kmeans/kmeans.c $(SRC_DIR)/kmeans/kmeans.h -o $@

$(OBJ_DIR)/kbag.o: $(SRC_DIR)/kbag/kbag.c 
	$(CC) $(CFLAGS) $(SRC_DIR)/kbag/kbag.c $(SRC_DIR)/kbag/kbag.h -o $@

$(OBJ_DIR)/gmmPval.o: $(SRC_DIR)/gmmPval/gmmPval.c 
	$(CC) $(CFLAGS) $(SRC_DIR)/gmmPval/gmmPval.c $(SRC_DIR)/gmmPval/gmmPval.h $(INC) -o $@

$(OBJ_DIR)/gaussMix.o: $(SRC_DIR)/gaussMix/gaussMix.c 
	$(CC) $(CFLAGS) $(SRC_DIR)/gaussMix/gaussMix.c $(SRC_DIR)/gaussMix/gaussMix.h $(INC) -o $@

$(OBJ_DIR)/gammaDist.o: $(SRC_DIR)/gammaDist/gammaDist.c 
	$(CC) $(CFLAGS) $(SRC_DIR)/gammaDist/gammaDist.c $(SRC_DIR)/gammaDist/gammaDist.h -o $@

$(OBJ_DIR)/classify.o: $(SRC_DIR)/classify/stat.c $(SRC_DIR)/classify/classify.c
	$(CC) $(CFLAGS) $(SRC_DIR)/classify/classify.h $(SRC_DIR)/classify/classify.c $(SRC_DIR)/classify/stat.c $(INC) -o $@

$(OBJ_DIR)/count.o: $(SRC_DIR)/chopKmers/count.c 
	$(CC) $(CFLAGS) $(SRC_DIR)/chopKmers/count.c $(SRC_DIR)/chopKmers/count.h  $(INC) -o $@

$(OBJ_DIR)/kmers.o: $(SRC_DIR)/chopKmers/kmers.c
	$(CC) $(CFLAGS) $(SRC_DIR)/chopKmers/kmers.c $(INC) -o $@

$(OBJ_DIR)/recordTasks.o: $(SRC_DIR)/chopKmers/recordTasks.c 
	$(CC) $(CFLAGS) $(SRC_DIR)/chopKmers/recordTasks.c $(INC) -o $@

$(OBJ_DIR)/minLength.o: $(SRC_DIR)/chopKmers/minLength.c
	$(CC) $(CFLAGS) $(SRC_DIR)/chopKmers/minLength.c $(INC) -o $@

$(OBJ_DIR)/interface.o: $(SRC_DIR)/interface/SAG_extrapolator.cpp
	$(CC) $(CFLAGS) $(SRC_DIR)/interface/SAG_extrapolator.cpp $(INC) -o $@

$(OBJ_DIR)/FastaParser.o: $(SRC_DIR)/interface/FastaParser.cpp $(SRC_DIR)/interface/Genome.cpp
	$(CC) $(CFLAGS) $(SRC_DIR)/interface/FastaParser.cpp $(INC) -o $@

sagex:
	g++ -Wall -g -o sagex $(SRC_DIR)/interface/SAG_extrapolator.cpp $(SRC_DIR)/interface/FastaParser.cpp $(SRC_DIR)/interface/Genome.cpp \
 $(SRC_DIR)/chopKmers/count.c $(SRC_DIR)/chopKmers/kmers.c $(SRC_DIR)/chopKmers/minLength.c $(SRC_DIR)/chopKmers/recordTasks.c \
 $(SRC_DIR)/gammaDist/gammaDist.c \
 $(SRC_DIR)/posixMat/matrix.c \
 $(SRC_DIR)/classify/stat.c \
 $(SRC_DIR)/gmmPval/gmmPval.c \
 $(SRC_DIR)/gaussMix/gaussMix.c \
 $(SRC_DIR)/kmeans/kmeans.c \
 $(SRC_DIR)/classify/classify.c \
 $(SRC_DIR)/kbag/kbag.c \
 -lpthread $(INC) 

clean:
	@rm -f $(OBJ_DIR)/*.o 
	@rm ./sagex

test: sagex
	@./sagex -i test/SAG.fasta -g test/metagenome.fasta
