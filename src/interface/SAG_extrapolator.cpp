#include "Fasta.hpp"
#include "helper.hpp"
#include "classify.h"
#include "matrix.h"
#include "gammaDist.h"
#include <getopt.h>

int main( int argc, char *argv[] ) {
    Options options;
    options.FetchArguments(argc, argv);
    options.ReviewArguments(argv);

    Fasta Metagenome(options.Gm);
    Metagenome.parse_fasta();
    
    if( options.seed > 0 )
	srand(options.seed) ;

    Fasta SAG(options.input);
    SAG.parse_fasta();

    int *hits = NULL; 
	int **hitPtr ; 
	if( options.kmerPCA != NULL )
		hitPtr = NULL ; 
	else
		hitPtr = &hits ;

    classify(&SAG, &Metagenome, &options, hitPtr);

    std::cout << "Finished successfully." << std::endl;
	return 0 ; 
}
