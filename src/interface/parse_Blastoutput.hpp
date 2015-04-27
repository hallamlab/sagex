#ifndef BLAST_H
#define BLAST_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <unistd.h>

using namespace std;

class BlastTable {
    ifstream blast_ptr;
    public:
        char **qseqid;
        char **sseqid;
        double *pident;
        int *length;
        int N_hits;
        BlastTable(char *); 
        int parseTable();
        int enterFields(char**, long int, long int, long int);
        int validateLine(int);
};

#endif
