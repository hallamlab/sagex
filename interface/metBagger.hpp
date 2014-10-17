#ifndef Gm_extract_H
#define Gm_extract_H

#include <cstdio>
#include <cstdlib>
#include <set>
#include <iostream>
#include <cstring>
#include "Genome.hpp"

using namespace std;

class metBagger : public Genome {
    public:
        int header_Nchar;
        metBagger();
        int record_header( string, int, int * );
        int record_sequence( string, int );
        int get_Gm_scaffolds(char **, double *, int *, int, char **, char **, int, int, int);
        int Bag( set<char *>, char **, char **, int );

        friend set<char *> find_unique_scaffolds(char ** , double *, int *, int, int, int);
};

#endif
