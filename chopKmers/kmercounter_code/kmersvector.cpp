#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include "utilities.h"
#include <vector>

 
int main( int argc, char **argv ){
    if( argc <= 1 ){
        std::cerr << "Usage: "<<argv[0]<<" [infile]" << std::endl;
        return -1;


    }
 
   std::string outputfile;
   std::vector<std::string> fasta_files;
   unsigned int size = 1000;

   for(int i = 1; i < argc ; i++) {
     //    std::cout << "Now after -c " << argv[i] << std::endl;
        if( strncmp(argv[i], "-f", strlen("-f")) == 0 ) {
           fasta_files.push_back(std::string(argv[++i]));
        }
        else if( strncmp(argv[i], "-o", strlen("-o")) == 0 ) {
           outputfile = argv[++i];
        }
        else if( strncmp(argv[i], "-s", strlen("-s")) == 0 ) {
           size = atol(argv[++i]);
        }
        else {
           std::cout << "ERROR: Choices for -f file name" << std::endl;
           return false;
        }
    }

    std::string line, name, content;
    int freq[256];
    
    std::ostream *output;
    std::ofstream fout;
    if( outputfile.size() > 0) {
        fout.open(outputfile.c_str(), std::ifstream::out);
        output = &fout;
    }


    for(std::vector<std::string>::iterator itfile = fasta_files.begin(); itfile != fasta_files.end(); itfile++) {
        std::ifstream input(itfile->c_str());
        if(!input.good()){
            std::cerr << "Error opening '"<<argv[1]<<"'. Bailing out." << std::endl;
            return -1;
        }
        while( std::getline( input, line ).good() ){
            if( line.empty() || line[0] == '>' ){ // Identifier marker
                if( !name.empty() ){ // Print out what we read from the last entry
          //          std::cout << name <<  content << std::endl;
                    computeKmerVectors(content, size, output);
                    name.clear();
                }
                if( !line.empty() ){
                    name = line.substr(1);
                }
                content.clear();
            } else if( !name.empty() ){
                if( line.find(' ') != std::string::npos ){ // Invalid sequence--no spaces allowed
                    name.clear();
                    content.clear();
                } else {
                    content += line;
                }
            }
        }
        if( !name.empty() ){ // Print out what we read from the last entry
         //   std::cout << name <<  content << std::endl;
            computeKmerVectors(content, size, output);
        }
        input.close();
    }
 
    return 0;
}
