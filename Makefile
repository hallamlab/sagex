
#===========================
# SAGex Makefile
# Connor Morgan-Lang 2015
# Hallam lab, Department of Microbiology and Immunology, UBC
#===========================

SHELL := /bin/bash -e

export OBJ_DIR = obj#For the compiled object files
export SRC_DIR = src#Where the c++ files live

CC=g++
CFLAGS=-Wall -g #removed '-c'
INC=-I$(SRC_DIR)/classify -I$(SRC_DIR)/posixMat -I$(SRC_DIR)/gammaDist -I$(SRC_DIR)/chopKmers -I$(SRC_DIR)/gaussMix -I$(SRC_DIR)/gmmPval -I$(SRC_DIR)/kbag -I$(SRC_DIR)/kmeans -I$(SRC_DIR)/interface
PROG=sagex

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
	$(OBJ_DIR)/Fasta.o \
	$(OBJ_DIR)/interface.o 

all: $(OBJ_DIR) $(PROG)

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

$(OBJ_DIR)/Fasta.o: $(SRC_DIR)/interface/Fasta.cpp
	$(CC) $(CFLAGS) $(SRC_DIR)/interface/Fasta.cpp $(INC) -o $@

$(PROG):
	@printf "\t\t Making sagex:\n"
	@printf "**************************************************\n"
	@printf "\n"
	$(CC) $(CFLAGS) -o $@ $(SRC_DIR)/interface/*.cpp \
 $(SRC_DIR)/chopKmers/*.c \
 $(SRC_DIR)/gammaDist/gammaDist.c \
 $(SRC_DIR)/posixMat/matrix.c \
 $(SRC_DIR)/classify/*.c \
 $(SRC_DIR)/gmmPval/gmmPval.c \
 $(SRC_DIR)/gaussMix/gaussMix.c \
 $(SRC_DIR)/kmeans/kmeans.c \
 $(SRC_DIR)/kbag/kbag.c \
 -lpthread $(INC) 

.PHONY: clean
clean:
	@rm -f $(OBJ_DIR)/*.o 
	@if [ -f $(PROG) ]; then rm $(PROG); else printf "$(PROG) has not been made.\n"; fi;

.PHONY: test
test: sagex
	./sagex -i test/mock_SAG.fasta -G test/mock_metaG.fasta \
-t 4 -v -K \
-X test/kmer20.test.tsv -o test/mock_output.test.fasta -Y test/kmerFrequencies.test.tsv
	./sagex -i test/mock_SAG.fasta -G test/mock_metaG.fasta \
-t 2 -v -k 8 -K \
-o test/mock_output.test.fasta
	./sagex -i test/mock_SAG.fasta -G test/mock_metaG.fasta -k 1 -o test/mock_output.test.fasta -v
	@rm test/kmer20.test.tsv test/mock_output.test.fasta test/kmerFrequencies.test.tsv
