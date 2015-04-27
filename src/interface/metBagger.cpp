#include "metBagger.hpp"
    
set<char *> find_unique_scaffolds(char ** sseqid, double *pidents, int *lengths, int N_hits, int pID, int minAL){
    set<char *> NR_scaffolds;
    NR_scaffolds.clear();
    int c;
    for (c = 0; c < N_hits; c++) {
        if ( (pidents[c] >= pID) && (lengths[c] >= minAL) )
            NR_scaffolds.insert(sseqid[c]);
        else
            ;
    }
    if ( NR_scaffolds.size() == 0 ) {
        cerr << "No significant hits were found in the Blast alignment!" << endl;
        cerr << "Please ensure the blast output table is not empty." << endl;
        exit(0);
    }
    return NR_scaffolds;
}

metBagger::metBagger() {
    header_Nchar = 0;
}

int metBagger::record_header(string S_header, int header_len, int *header_Nchar) {
    * header_Nchar += header_len;
    header = (char **) realloc (header, * header_Nchar * sizeof(char *));
    header[N_contigs] = (char *) malloc (header_len * sizeof(char));
    strcpy(header[N_contigs], S_header.c_str());
    return 0;
}

int metBagger::record_sequence(string S_seq, int contig_len) {
    sequence = (char **) realloc (sequence, (genome_length + 1) * sizeof(char *));
    sequence[N_contigs] = (char *) malloc ( (genome_length + (1*N_contigs)) * sizeof(char));
    strcpy(sequence[N_contigs], S_seq.c_str());
    N_contigs++;
    return 0;
}

int metBagger::Bag(set<char *> NR_scaffolds, char **MetaG_headers, char **MetaG_sequences, int MetaG_N_contigs){
    int n;
    int header_len, contig_len;
    string S_header, S_seq;
    if (NR_scaffolds.empty() == true) {
        for(n=0; n < MetaG_N_contigs; n++) {
if( n % 10000 == 0 ){ printf( "DEBUG n: %i, goal: %i\n" , n , MetaG_N_contigs ) ; } 
            header_len = strlen(MetaG_headers[n]) + 2;
            S_header = string(MetaG_headers[n]);
            record_header(S_header, header_len, &header_Nchar);
            
            contig_len = strlen(MetaG_sequences[n]);
            genome_length += contig_len;
            S_seq = string(MetaG_sequences[n]);
            record_sequence(S_seq, contig_len);
        }
    }
    else {
        char *blast_header = (char *) malloc (2 * sizeof(char));
        strcpy(blast_header, ">");
        for (set<char *>::iterator it=NR_scaffolds.begin(); it!=NR_scaffolds.end(); ++it) { 
            blast_header = (char *) realloc (blast_header, (strlen(*it) + 1) * sizeof(char));
            strcat(blast_header, *it);
            for(n=0; n < MetaG_N_contigs; n++) {
                if (strcmp( blast_header, MetaG_headers[n]) == 0) {
                    header_len = strlen(MetaG_headers[n]) + 2;
                    S_header = string(MetaG_headers[n]);
                    record_header(S_header, header_len, &header_Nchar);
        
                    contig_len = strlen(MetaG_sequences[n]);
                    genome_length += contig_len;
                    S_seq = string(MetaG_sequences[n]);
                    record_sequence(S_seq, contig_len);
                }
                else {
                    ;
                }
            }
            strcpy(blast_header, ">");
        }
        if (N_contigs == 0 ) {
            cerr << "No metagenomic scaffolds in the blast output table were found in the metagenome assembly selected." << endl;
            cerr << "Please make sure the blast table file  corresponds to the metagenome provided!" << endl;
            exit(0);
        }
    }
    return 0;
}

int metBagger::get_Gm_scaffolds(char **sseqid, double *pidents, int *lengths, int N_hits, char **headers, char **sequences, int MetaG_N_contigs, int pID, int minAL) {
    set<char *> NR_scaffolds;
    if (sseqid != NULL) {
        NR_scaffolds = find_unique_scaffolds(sseqid, pidents, lengths, N_hits, pID, minAL);
        free(sseqid);
    }
    Bag(NR_scaffolds, headers, sequences, MetaG_N_contigs);
    exit(0) ; 
}
