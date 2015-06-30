#include "helper.hpp"

void Options::Print_help() {
    std::string help;
    help = "\nHelp description:\n\
    -i: \n \
    \tThe SAG assembly. This must be in FASTA format \n\
    -G: \n \
    \tThe Metagenome assembly to be used for genome extrapolation. This must be in FASTA format \n\
    -o \n\
    \tThe output file. If not specified, output sequences in FASTA format are printed to stdout [OPTIONAL]\n\
    \n\t\tContig uniformity options:\n\
    -c: \n\
    \tThe desired length (in base-pairs) of the chopped contig. If less than one, entire contigs will be used. [DEFAULT is proportions]\n\
    -x: \n\
    \tThe number of base-pairs overlapping between each of the chopped contigs. [DEFAULT is proportions]\n\
    \n\t\tClassification and Extrapolation-specific options:\n\
    -k: \n \
    \tDesired number of Gaussians in the mixture model. [DEFAULT = 2]\n\
    -K: \n \
    \tFlag toggling whether to reduce number of Gaussians if needed [DEFAULT = False]\n\
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
    -s: \n\
    \tThe seed for random number generation. [DEFAULT = -1]\n\
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
    std::cerr << help << std::endl;
    return;
}

void Options::Print_usage(char* arg) {
    std::cerr << "\nUSAGE: \t" << arg << " [options] -i sag.fasa -G metagenome.fasta" << std::endl;
    return;
}

void Options::ReviewArguments(char* argv[]) {
    if (this->hFlag == true) {
        this->Print_usage(argv[0]);
        this->Print_help();
        exit(2);
    }
    if (this->problem == true) {
        exit(3);
    }
    if ( (!proportionFlag && overlap == 0) || (!proportionFlag && chopSize) ){
        std::cerr << "Only chopSize or overlap were specified when both are required." << std::endl;
        this->Print_help();
        exit(4);
    }
    else if (errFlag > 0){
        std::cerr << "Mutually exclusive arguments were included in the command-line options." << std::endl;
        this->Print_help();
        exit(5);
    }
    else if ((!input) || (!Gm)) {
        std::cerr << "Required arguments";
        if (input == NULL)
            std::cerr << " -i ";
        if (Gm == NULL)
            std::cerr << " -G ";
        std::cerr << "have not been included!" << std::endl;
        this->Print_usage(argv[0]);
        exit(6);
    }
    return;
}

int Options::FetchArguments(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Type '" << argv[0] << " -h' for help options." << std::endl;
        this->problem = true;
        return 1;
    }
    int c;
    unsigned int j;
    char const * arguments_need_input = "iGocxABEmMtkXYCDs";
    while ((c = getopt (argc, argv, "i:G:o:c:x:A:B:E:m:M:t:k:X:Y:C:D:s:hvK")) != -1) {
        switch (c) {
            case 'i':
                this->input = optarg;
                break;
            case 'G':
                this->Gm = optarg;
                break;
            case 'D':
                this->mahalaMultiple = atof(optarg) ;
                break;
            case 'K':
                this->fixK = true;
                break;
            case 'C':
                this->lcsCut = atoi(optarg) ;
                break;
            case 'o':
                this->output = optarg;
                break;
            case 's':
                this->seed = atoi(optarg) ;
                break ;
            case 'c':
                this->chopSize = atoi(optarg);
                this->proportionFlag = false;
                break;
            case 'x':
                this->overlap = atoi(optarg);
                this->proportionFlag = false;
                break;
            case 'A':
                this->Alpha = atof(optarg);
                break;
            case 'B':
                this->Beta = atof(optarg);
                break;
            case 'E':
                this->eps = atof(optarg);
                break;
            case 'm':
                this->itLen = atoi(optarg);
                break;
            case 'M':
                this->itMax = atoi(optarg);
                break;
            case 't':
                this->threads = atoi(optarg);
                break;
            case 'k':
                this->k = atoi(optarg);
                break;
            case 'X':
                this->kmerPCA = optarg;
                break;
            case 'Y':
                this->kmerFreq = optarg;
                break ;
            case 'v':
                this->verbose = true;
                break;
            case 'h':
                this->hFlag = true;
                break;
            case '?':
                for (j=0; j < strlen(arguments_need_input); j++) {
                    if (optopt == arguments_need_input[j]) {
                        fprintf(stderr, "Option -%c requires an argument. \n", optopt);
                        exit(1);
                    }
                }
                if (isprint (optopt)) {
                    fprintf(stderr, "Unknown option character -%c. \n", optopt);
                    this->problem = true;
                }
                else {
                    fprintf(stderr, "Unknown option `\\x%x'.\n", optopt);
                    this->problem = true;
                }
                return 1;
            default:
                exit(0);
        }
    }
    return 0;
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}