#!/bin/bash

cd $1

./../sagex -i sag.fasta -G mock_metaG.fasta -t 4 -v -P -c -1 -C 20 -X kmer20.tmp -o fasta20.tmp -Y kmerFrequencies_20.tsv
