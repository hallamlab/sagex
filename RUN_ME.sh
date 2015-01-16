#!/bin/bash
./sagex -i ../data/test_data/sag.fasta -G ../data/test_data/metaG.fasta -v -P -c -1 -C 20 -X ./tmp/kmer20 > ./tmp/fasta20
