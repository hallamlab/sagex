sagex
=====

$\Sigma$

##Purpose:
SAGex is a lightweight genome extrapolation tool. It will take any genomic nucleotide sequence (in FASTA format) and a second nucleotide sequence, and recruit all contigs that are in the first
genome's tetranucleotide space. This was originally intended for Single-cell Amplified Genomes (SAGs) and Metegenomes for recruitment. 

##Installation:
Decompress the compressed sagex repository, if needed. 
```
cd SAGex
make
```
This will compile the source files and create the executable `sagex`. To make sure everything has compiled properly on your machine try:
```
make test
```

Refer to the [Issues]() page for known errors during installation and run-time.

##Citation:
Please use:
Durno, W.E., Morgan-Lang, C., Hawley, A.K., and Hallam, S.J. SAGex. Bioinformatics. 201X.
