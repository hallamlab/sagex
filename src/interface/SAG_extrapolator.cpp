//#include "parse_Blastoutput.hpp"
//#include "metBagger.hpp"
#include "FastaParser.hpp"
#include "Genome.hpp"
#include "classify.h"
#include "matrix.h"
#include "gammaDist.h"
#include <getopt.h>

void checkArgs(int hFlag, int errFlag, char *input, char *Gm, char *blastout, int argc) {
    string usage = "\nUSAGE: \t ./sagex [options] -i sag.fasta -G metagenome.fasta";
    string help = "\nHelp description:\n\
-i: \n \
\tThe SAG assembly. This must be in FASTA format \n\
-o \n\
\tThe output file. If not specified, output sequences in FASTA format are printed to stdout [OPTIONAL]\n\
-G: \n \
\tThe Metagenome assembly to be used for genome extrapolation. This must be in FASTA format \n\
-b: \n \
\tThe BLAST alignment file of the SAG to the Metagenome. This must be in output format 6 [OPTIONAL]\n\
\n\t\tBLAST-specific options:\n\
-p: \n \
\tThe minimum percentage of identical sequence between the SAG and metagenomic contigs for the contig to be included in the analysis. Integer type. [DEFAULT = 85]\n\
-a: \n \
\tAn integer representing the minimum number of base-pairs of the SAG aligned to the metagenome to be considered a good hit. [DEFAULT = 2000]\n\
\n\t\tContig uniformity options:\n\
-c: \n\
\tThe desired length (in base-pairs) of the chopped contig. If less than one, entire contigs will be used. [DEFAULT = 2000]\n\
-x: \n\
\tThe number of base-pairs overlapping between each of the chopped contigs. [DEFAULT = 500]\n\
-P: \n\
\tDo not chop contigs. Run kmer proportions. Supersedes -x.\n\
\n\t\tClassification and Extrapolation-specific options:\n\
-k: \n \
\tDesired number of Gaussians in the mixture model. [DEFAULT = 1]\n\
-A: \n\
\tThe false positive rate for attributing contig chops to SAGs. [DEFAULT = 0.05]\n\
-B: \n\
\tThe proportion of contig chops required before deciding the contig is in the SAG. [DEFAULT = 0.8]\n\
-E: \n\
\tConvergance parameter for Eigen value calculation. [DEFAULT = 0.0000001]\n\
-m: \n\
\tMinimum contig length. [DEFAULT = 2000]\n\
-M: \n\
\tThe maximum number of iterations for numerical integration and bootstraps. [DEFAULT = 1000]\n\
\n\t\tOther options:\n\
-X: \n\
\tA file to write the kmer PCA in an R-readable format.\n\
-Y: \n\
\tA file to report the kmer frequencies in an R-readable format.\n\
-v: \n\
\tTurns on verbosity.\n\
-t \n\
\tThe maximum number of POSIX threads allowed. [DEFAULT = 1]\n\
-h \n\
\tShow help message and exit\n";
    if (hFlag == 1) {
        cerr << usage << endl << help << endl;
        exit(1);
     }
    else if (errFlag > 0){
        string ME_error = "Mutually exclusive arguments were included in the command-line options.";
        cerr << ME_error << endl << help << endl;
        exit(2);
    }
    else if ((!input) || (!Gm)) {
        cerr << "Required arguments";
        if (input == NULL)
            cerr << " -i ";
        if (Gm == NULL)
            cerr << " -G ";
        cerr << "have not been included!" << endl << usage << endl;
        exit(3);
    }
    //exit(0) ;  
}   
/*
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
*/
int main( int argc, char *argv[] ) {
    
    char *input = NULL;
    char *Gm = NULL;
    char *blastout = NULL;
    char *output = NULL; 
    char *kmerPCA = NULL; 
    char *kmerFreq = NULL; 
    int k = 1;
    //int pID = 85;
    //int minAL = 2000;
    int chopSize = 2000;
    int overlap = 500;
    double Alpha = 0.05;
    double Beta = 0.8;
    double eps = 0.0000001; 
    double D = 1.0 ; // malahanobis multiple  
    int itLen = 2000;
    int itMax = 10000; 
    int lcsCut = -1 ; 
    int threads = 1;
    int seed = -1 ; 
    int K = -1 ; 
    int c, errFlag, verbose, proportionFlag, hFlag;
    errFlag = verbose = proportionFlag = hFlag = opterr = 0;
    while ((c = getopt (argc, argv, "i:G:b:o:c:x:A:B:E:m:M:t:k:X:Y:C:D:s:hPvK")) != -1) {
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
	    case 'D': 
		D = atof(optarg) ; 
		break; 
	    case 'K':
		K = 1; 
		break; 
/*
            case 'p':
                pID = atoi(optarg);
                break;
            case 'a':
                minAL = atoi(optarg);
                break;
*/
            case 'C':
                lcsCut = atoi(optarg) ; 
                break;
            case 'o':
                output = optarg;
                break;
	    case 's':
	        seed = atoi(optarg) ; 
		break ; 
            case 'c':
                    chopSize = atoi(optarg);
                break;
            case 'x':
                if (proportionFlag == 1) {
                    proportionFlag = -1;
                    errFlag++;
                }
                else
		            overlap = atoi(optarg);
		        break;
            case 'P':
                if (proportionFlag != 0)
                    errFlag++;
                else
                    proportionFlag++;
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
                itLen = atoi(optarg);
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
                kmerPCA = optarg;
                break;
            case 'Y':
                kmerFreq = optarg; 
                break ; 
            case 'v':
                verbose = 1;
                break;
            case 'h':
                hFlag = 1;
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
    
    checkArgs(hFlag, errFlag, input, Gm, blastout, argc);
	
    FastaParser Metagenome(Gm); 
    Metagenome.parse_fasta(); 
    //cout << "Metagenome size = " << Metagenome.genome_length << endl; 
    
    if( seed > 0 ) 
	srand(seed) ; 

    FastaParser SAG(input); 
    SAG.parse_fasta(); 
    //cout << "SAG size = " << SAG.genome_length << endl;
/*    
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
	if( kmerPCA != NULL ) 
		hitPtr = NULL ; 
	else
		hitPtr = &hits ;

    classify ( SAG.sequence , SAG.N_contigs , SAG.header , Metagenome.sequence , Metagenome.N_contigs , Metagenome.header , Alpha , Beta , threads , eps , itLen , itMax , chopSize , overlap , proportionFlag , k , K , D , verbose , kmerFreq , kmerPCA , hitPtr , output , lcsCut );
	 
	return 0 ; 
}
