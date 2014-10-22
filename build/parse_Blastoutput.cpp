#include "parse_Blastoutput.hpp"

BlastTable::BlastTable( char * input ) {
    blast_ptr.open(input);
    if ( !blast_ptr.is_open()) {
        fprintf(stderr, "Could not open file %s.\n Exiting now!\n", input);
        exit(0);
    }
    else {
        qseqid = (char **) malloc (1 * sizeof(char*));
        sseqid = (char **) malloc (1 * sizeof(char*));
        pident = (double *) malloc (1 * sizeof(double));
        pident[0] = 1;
        length = (int *) malloc (1 * sizeof(int));
        length[0] = 1;
    }
}

char **str_slicer(char *s, const char *delim, int * slices_size = 0) {
    string t = string(s);
    string substring;
    int delim_size = strlen(delim);
    int str_len = t.length();
    char **slices = (char **) malloc ((str_len + 1) * sizeof(char*));
    int x, n, cur_len;
    n = 0;
    slices[n] = (char *) malloc (1 * sizeof(char));
    slices[n] = NULL;
    for (x = 0; x < str_len; x++) {
        substring = t.substr(x, delim_size);
        if (strcmp(substring.c_str(), delim) == 0) {
            n++;
            slices[n] = (char *) malloc (1 * sizeof(char));
            slices[n] = NULL;
        }
        else {
            if ( slices[n] != NULL )
                cur_len = strlen(slices[n]);
            else
                cur_len = 0;
            slices[n] = (char *) realloc (slices[n], (cur_len + 2) * sizeof(char));
            slices[n][cur_len] = t.at(x);
            slices[n][cur_len + 1] = '\0';
        }
    }
    *slices_size = n+1;
    return slices;
}

int BlastTable::enterFields(char **fields, long int x, long int query_Nchar, long int subject_Nchar) {
        qseqid = (char **) realloc (qseqid, query_Nchar * sizeof(char*));
//        cerr << "DEBUG 1" << endl;
        qseqid[x] = (char *) malloc ((strlen(fields[0]) + 1) * sizeof(char));
        strcpy(qseqid[x], fields[0]);
        sseqid = (char **) realloc (sseqid, subject_Nchar * sizeof(char*));
//        cerr << "DEBUG 2" << endl;
        sseqid[x] = (char *) malloc ((strlen(fields[1]) + 1) * sizeof(char));
        strcpy(sseqid[x], fields[1]);
//        cerr << "DEBUG 3" << endl;
        pident[0]++;
        pident = (double *) realloc (pident, (pident[0] + 1) * sizeof(double));
        pident[x+1] = atof(fields[2]);
//        cerr << "DEBUG 4" << endl;
        length[0]++;
        length = (int *) realloc (length, (length[0] + 1) * sizeof(int));
        length[x+1] = atoi(fields[3]);
//        cerr << "DEBUG 5" << endl;
        return 0;
}

int BlastTable::parseTable() {
    string s;
    const char tab[2] = "\t";
    char **fields;
    char * line = (char *) malloc (1 * sizeof(char));
    long int x, str_len, query_Nchar, subject_Nchar;    
    int line_Nelements = 0; 
    if (blast_ptr.good())
        getline( blast_ptr, s);
    else { 
        fprintf(stderr, "Something is wrong with the blast output table. Exiting now!\n");
        exit(0);
    }
    x = 0;
    query_Nchar = subject_Nchar = 0;
    while( blast_ptr.good() ) {
        str_len = s.length();
        line = (char *) realloc (line, str_len + 1 * sizeof(char));
        strcpy(line, s.c_str());
        fields = str_slicer(line, tab, &line_Nelements);
        query_Nchar += (strlen(fields[0]) + 2);
        subject_Nchar += (strlen(fields[1]) + 2);
        validateLine(line_Nelements);
        enterFields(fields, x, query_Nchar, subject_Nchar);
        x++;
        getline( blast_ptr, s);
    }
    free(fields);
    N_hits = x;
}

int BlastTable::validateLine(int line_Nelements) {
    string warning;
    warning = "The blast output table was not formatted correctly!\n\
\tThere should be 12 fields by using the blastn command outfmt 6.\n\
\tPlease reformat your blast alignment table and try again.\n";
    if ( line_Nelements == 12 ) {
        return 0;
    }
    else {
        fprintf(stderr, "A line contains only %d fields.\n", line_Nelements);
        cerr << warning << endl;
        exit(0);
    }
}
