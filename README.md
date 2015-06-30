sagex
=====

##Purpose:
SAGex is a lightweight genome extrapolation tool. It will take any genomic nucleotide sequence (in FASTA format) and a second nucleotide sequence, and recruit all contigs that are in the first
genome's tetranucleotide space. This was originally intended for Single-cell Amplified Genomes (SAGs) and Metegenomes for recruitment. 

##Installation:
Decompress the compressed sagex repository, if needed. 
```
git clone git@github.com:hallamlab/sagex.git
cd sagex
make sagex
```
This will compile the source files and create the executable `sagex`. 
To make sure everything has compiled properly on your machine try:
```
make test
```
There should be no errors with this test.

Refer to the [Issues](https://github.com/hallamlab/sagex/issues) page for known errors during
installation and run-time. If you have encountered a new error please make an issue with the
appropriate label.

##Citation:
This is an academic software so when you use sagex please cite using the following:
Durno, W.E., Hawley, A.K., Morgan-Lang, C., and Hallam, S.J. SAGex. Bioinformatics. 201X.
