#include "parse_Blastoutput.hpp"
#include "FastaParser.hpp"
#include "metBagger.hpp"
#include "Genome.hpp"
#include "classify.h"
#include "matrix.h"
#include "gammaDist.h"

int checkArgs(int hflag, int argc) {
    string usage = "USAGE: \t ./extrapolate_SAG [options] -i sag.fasta -G metagenome.fasta -b table.blastout\n\n\
-h\t\tShow help message and exit\n\
-o\t\tThe output Directory [OPTIONAL]\n"; /*Needs editing*/
    if (hflag == 1) {
        string help = "\nHelp description:\n\
-i: \n \
\tThe SAG assembly. This must be in FASTA format \n\
-G: \n \
\tThe Metagenome assembly to be used for genome extrapolation. This must be in FASTA format \n\
-b: \n \
\tThe BLAST alignment file of the SAG to the Metagenome. This must be in output format 6 \n\
\n\t\tBLAST-specific options:\n\
-p: \n \
\tThe minimum percentage of identical sequence between the SAG and metagenomic contigs for the contig to be included in the analysis. [DEFAULT = 85]\n\
-a: \n \
\tAn integer representing the minimum number of base-pairs of the SAG aligned to the metagenome to be considered a good hit. [DEFAULT = 2000]\n\
-k: \n \
\tDesired number of Gaussians in the mixture model. [DEFAULT = 1]\n\
\n\t\tContig uniformity options:\n\
-c: \n\
\tThe desired length (in base-pairs) of the chopped contig. [DEFAULT = 2000]\n\
-x: \n\
\tThe number of base-pairs overlapping between each of the chopped contigs. [DEFAULT = 500]\n\
\n\t\tClassification and Extrapolation-specific options:\n\
-A: \n\
\tThe false positive rate for attributing contig chops to SAGs. [DEFAULT = 0.05]\n\
-B: \n\
\tThe proportion of contig chops required before deciding the contig is in the SAG. [DEFAULT = 0.8]\n\
-E: \n\
\tConvergance parameter for Eigen value calculation. [DEFAULT = 0.0001]\n\
-m: \n\
\tThe minimum number of iterations for Eigen value calculation. [DEFAULT = 10]\n\
-M: \n\
\tThe maximum number of iterations for Eigen value calculation. [DEFAULT = 1000]\n\
\n\t\tOther options:\n\
-X: \n\
\tReport kmers in an R-readable format instead of hit names.\n\
-v: \n\
\tTurns on verbosity.\n\
-t \n\
\tThe maximum number of POSIX threads allowed. [DEFAULT = 1]\n";
        cerr << usage << endl << help << endl;
        exit(1);
     }
    else if (argc < 4) {
        fprintf(stderr, "Required arguments have not been included\n\n");
        cerr << usage << endl;
        exit(0);
    }
}   

int *map_chops(char ** metBag_header, int metBag_Ncontigs) {
    int *metBag_headerMap = (int *) malloc (metBag_Ncontigs * sizeof(int));
    int c, n;
    n = 0;
    char *lastHeader = metBag_header[0];
    for (c = 0; c < metBag_Ncontigs; c++) {
        if (strcmp(lastHeader, metBag_header[c]) != 0) {
            strcpy(lastHeader, metBag_header[c]);
            n++;
        }
        else
            ;
        metBag_headerMap[c] = n;
    }
    return metBag_headerMap;
}

int main( int argc, char *argv[] ) {
    
    char *input = NULL;
    char *Gm = NULL;
    char *blastout = NULL;
    char *output = NULL; 
    int k = 1; 
    int pID = 85;
    int minAL = 2000;
    int hflag = 0;
    int chopSize = 2000;
    int overlap = 500;
    double Alpha = 0.05;
    double Beta = 0.8;
    double eps = 0.0001;
    int itMin = 10;
    int itMax = 100000;
    int threads = 1;
    int kmerFlag = 0 ; 
    int verbose = 0 ; 
    int c;
    opterr = 0;
    while ((c = getopt (argc, argv, "i:G:b:o:p:a:c:x:A:B:E:m:M:t:k:hXv")) != -1) {
        switch (c) {
            case 'i':
                input = optarg;
                break;
            case 'G':
                Gm = optarg;
                break;
            case 'b':
                blastout = optarg;
                break;
            case 'p':
                pID = atoi(optarg);
                break;
            case 'a':
                minAL = atoi(optarg);
                break;
            case 'o':
                output = optarg;
                break;
            case 'c':
                chopSize = atoi(optarg);
                break;
            case 'x':
                overlap = atoi(optarg);
                break;
            case 'A':
                Alpha = atof(optarg);
                break;
            case 'B':
                Beta = atof(optarg);
                break;
            case 'E':
                eps = atof(optarg);
                break;
            case 'm':
                itMin = atoi(optarg);
                break;
            case 'M':
                itMax = atoi(optarg);
                break;
            case 't':
                threads = atoi(optarg);
                break;
            case 'k':
                k = atoi(optarg);
                break;
            case 'X':
                kmerFlag = 1;
                break;
            case 'v':
                verbose = 1;
                break;
            case 'h':
                hflag = 1;
                break;
            case '?':
                if ((optopt == 'i') || (optopt == 'G') || (optopt == 'o') || (optopt == 'c') || (optopt == 'x') || (optopt == 'A') || (optopt == 'B') || (optopt == 't')) {
                    fprintf(stderr, "Option -%c requires an argument. \n", optopt);
                    exit(1);
                }
                else if (isprint (optopt))
                    fprintf(stderr, "Unknown option character -%c. \n", optopt);
                else
                    fprintf(stderr, "Unknown option `\\x%x'.\n", optopt);
                return 1;
            default:
                exit(0);
        }
    }
    
    checkArgs(hflag, argc);
	
    FastaParser Metagenome(Gm); 
    Metagenome.parse_fasta(); 
    //cout << "Metagenome size = " << Metagenome.genome_length << endl;

    FastaParser SAG(input); 
    SAG.parse_fasta(); 
    //cout << "SAG size = " << SAG.genome_length << endl;
    
    metBagger metBag;
    
    if (blastout != NULL) {

        BlastTable blastclass(blastout); 
        blastclass.parseTable(); 
    
        free(blastclass.qseqid); 

        metBag.get_Gm_scaffolds(blastclass.sseqid, blastclass.pident, blastclass.length, blastclass.N_hits, Metagenome.header, Metagenome.sequence, Metagenome.N_contigs, pID, minAL);
	
        //cout << "MetBag size = " << metBag.genome_length << endl;
    }
    else {
        // metBag.get_Gm_scaffolds(NULL, NULL, NULL, 0, Metagenome.header, Metagenome.sequence, Metagenome.N_contigs, pID, minAL);
        metBag.header = Metagenome.header ; 
	metBag.sequence = Metagenome.sequence ; 
	metBag.N_contigs = Metagenome.N_contigs ; 
	metBag.genome_length = Metagenome.genome_length ; 
	// cout << "Passed get_Gm_scaffolds" << endl;
    }
	/*
printf( "DEBUG chopping Gm...\n" ) ; 
        Chopper chopped_metBag(chopSize, overlap);
        chopped_metBag.chop_scaffolds(metBag.header, metBag.sequence, metBag.N_contigs);

printf( "DEBUG chopping SAG...\n" ) ; 
        Chopper chopped_SAG(chopSize, overlap);
        chopped_SAG.chop_scaffolds(SAG.header, SAG.sequence, SAG.N_contigs);
	
        int *metBag_headerMap = map_chops(chopped_metBag.header, chopped_metBag.N_contigs);
    	*/
	
        int *hits = NULL; 
	int **hitPtr ; 
	if( kmerFlag > 0 ) 
		hitPtr = NULL ; 
	else
		hitPtr = &hits ; 
        classify ( SAG.sequence , SAG.N_contigs , metBag.sequence , metBag.N_contigs , Alpha , Beta , threads , eps , itMin , itMax , chopSize , overlap , k , verbose , hitPtr );
	
	if( kmerFlag > 0 ) 
		return 0 ; 
	
	for( int i = 0 ; i < metBag.N_contigs ; i++ ) 
	{
		if( hits[i] > 0 ) 
		{
			cout << metBag.header[i] << endl ; 
			cout << metBag.sequence[i] << endl ; 
		}
	} 
	return 0 ; 
}
